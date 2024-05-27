#!/usr/bin/python3.10

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
    num_frames = len(u.trajectory)   #[args.bf:args.ef:args.dtf])
    num_atoms = max([int(atom.name[2:]) for atom in u.select_atoms(selection).split('residue')[0]])

    opu_list = [[] for _ in range(num_atoms - 1)]
    with tqdm(total=num_frames, desc="Processing Frames") as pbar:
        for ts in u.trajectory:     #[args.bf:args.ef:args.dtf]:
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
                            opu_list[n].append(order_parameter)

            pbar.update(1)

    opu_array = list_to_array(opu_list)
    return abs(np.nanmean(opu_array, axis=1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculation of order parameters')
    parser.add_argument('-s', type=str, help='Topology file (default: %(default)s)', default='topol.tpr',
                        metavar='<.tpr>')
    parser.add_argument('-f', type=str, help='Trajectory file (default: %(default)s)', default='traj_comp.xtc',
                        metavar='<.xtc>')
    parser.add_argument('-o', type=str, help='Output text file (default: %(default)s)',
                        default='order_params.txt', metavar='<.dat/.txt/...>')
    parser.add_argument('-l', '--lipids', type=str, help='Lipid model used',
                        choices=['popc', 'dopc', 'dmpc'], required=True)
    parser.add_argument('-d', '--distance', type=float, metavar='<distance>', default=5,
                        help='Distance for determining the nearest lipids (default: %(default)sA)')
    parser.add_argument('--nearest', help='Turn on calculation of nearest lipids', action='store_true')
    parser.add_argument('--refmol', type=str,
                        help='Reference molecule around which the environment is studied (default: %(default)s)',
                        default='MOL', metavar='<resname>')
    # parser.add_argument('-bf', type=int, metavar='<start frame>', default=0,
    #                     help='Frame to start analysis from (default: %(default) frame)')
    # parser.add_argument('-ef', type=int, metavar='<end frame>', default=-1,
    #                     help='End frame for analysis (default: %(default) frame)')
    # parser.add_argument('-dtf', type=float, metavar='<frame step>', default=1,
    #                     help='Frame step to calculate order parameter (default: %(default))')
    # parser.add_argument('--selection', type=str, help='Another atom selection, can be used if lipid name')

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
