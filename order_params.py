#!/usr/bin/env python3

import numpy as np
import MDAnalysis as mda
import numpy.linalg as norm
from tqdm import tqdm
import argparse


def list_to_array(sequence):
    max_length = max(len(s) for s in sequence)
    padded_sequence = np.array(
        [np.pad(s, (0, max_length - len(s)), 'constant', constant_values=np.nan) for s in sequence])
    return padded_sequence


def order_parameters(selection):
    time = u.trajectory.totaltime

    if args.b is not None or args.e is not None or args.ts is not None:
        total_frames = len(u.trajectory)
        frames_range = range(int(args.b * total_frames / time) if args.b is not None else 0,
                             int(args.e * total_frames / time) if args.e is not None else len(u.trajectory) - 1,
                             int(args.ts * total_frames / time) if args.ts is not None else 1)
        frames = [i for i in frames_range]
    else:
        frames_range = range(args.bf if args.bf is not None else 0,
                             args.ef if args.ef is not None else len(u.trajectory),
                             args.step if args.step is not None else 1)
        frames = [i for i in frames_range]

    num_frames = len(u.trajectory[frames])

    num_atoms = max([int(atom.name[2:]) for atom in u.select_atoms(selection).split('residue')[0]])

    opu_list = [[] for _ in range(num_atoms - 1)]
    with tqdm(total=num_frames, desc="Processing Frames") as pbar:
        for ts in u.trajectory[frames]:
            lipid_atoms = u.select_atoms(selection).split('residue')
            for lipid in lipid_atoms:
                for atom in lipid.atoms:
                    n = int(atom.name[2:]) - 2
                    for h in atom.bonded_atoms:
                        if 'H' in h.name:
                            # Calculate C-H vector 
                            vector = h.position - atom.position

                            # Calculate cosine theta
                            cos_theta = np.dot(vector, [0, 0, 1]) / (np.linalg.norm(vector) * np.linalg.norm([0, 0, 1]))

                            # Calculate order parameter and append it to list
                            order_parameter = 0.5 * (3 * cos_theta ** 2 - 1)
                            opu_list[n].append(round(order_parameter, 3))

            pbar.update(1)

    opu_array = list_to_array(opu_list)
    return abs(np.nanmean(opu_array, axis=1))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculation of order parameters')

    parser.add_argument('-s', type=str,
                        help='Input topology file used to define the system structure (default: %(default)s)',
                        default='topol.tpr',
                        metavar='<.tpr>')

    parser.add_argument('-f', type=str,
                        help='Input trajectory file containing the molecular dynamics simulation data (default: %('
                             'default)s)',
                        default='traj_comp.xtc',
                        metavar='<.xtc>')

    parser.add_argument('-o', type=str,
                        help='Output file for writing calculated order parameters (default: %(default)s)',
                        default='order_params.txt',
                        metavar='<.dat/.txt/...>')

    parser.add_argument('-l', '--lipids', type=str,
                        help='Lipid model used in the simulation (required)',
                        choices=['popc', 'dopc', 'dmpc'],
                        required=True)

    parser.add_argument('-d', '--distance', type=float,
                        metavar='<distance>',
                        default=5.0,
                        help='Threshold distance in Angstroms for identifying nearest lipids (default: %(default)s Å)')

    parser.add_argument('--nearest',
                        help='Enable calculation of nearest lipids',
                        action='store_true')

    parser.add_argument('--refmol', type=str,
                        help='Residue name of the reference molecule around which the environment is studied ('
                             'default: %(default)s)',
                        default='MOL',
                        metavar='<resname>')

    parser.add_argument('-bf', type=int,
                        metavar='<first frame>',
                        default=None,
                        help='Frame number to start the analysis from (default: 0)')

    parser.add_argument('-b', type=int,
                        metavar='<start time>',
                        default=None,
                        help='Start time in picoseconds (ps) to begin the analysis (default: 0 ps)')

    parser.add_argument('-ef', type=int,
                        metavar='<last frame>',
                        default=None,
                        help='Frame number to end the analysis at')

    parser.add_argument('-e', type=int,
                        metavar='<last time>',
                        default=None,
                        help='End time in picoseconds (ps) to stop the analysis')

    parser.add_argument('-step', type=int,
                        metavar='<step>',
                        default=None,
                        help='Interval between frames to analyze (default: 1)')

    parser.add_argument('-ts', type=int,
                        metavar='<time step>',
                        default=None,
                        help='Time step interval in picoseconds (ps) for the analysis')

    args = parser.parse_args()

    u = mda.Universe(args.s, args.f)

    if args.nearest:
        sat = f'(resname {args.lipids.upper()} and name C3* and not name C3) and same resid as around' \
              f' {args.distance} resname {args.refmol}'
        unsat = f'(resname {args.lipids.upper()} and name C2* and not name C2) and same resid as around' \
                f' {args.distance} resname {args.refmol}'
    else:
        sat = f'resname {args.lipids.upper()} and name C3* and not name C3'
        unsat = f'resname {args.lipids.upper()} and name C2* and not name C2'

    opu = order_parameters(unsat)
    ops = order_parameters(sat)

    dif = len(opu) - len(ops)
    ops_extended = np.pad(ops, (0, dif), 'constant', constant_values=np.nan)
    result = np.column_stack((ops_extended, opu))
    np.savetxt(args.o, result, header='Saturated Unsaturated', comments='')
