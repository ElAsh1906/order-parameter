#!/usr/bin/env python3

import numpy as np
import MDAnalysis as mda
import numpy.linalg as norm
from tqdm import tqdm
import argparse
import pickle


def select_lipids(args):
    sat = f'resname {args.lipids.upper()} and name C3* and not name C3'
    unsat = f'resname {args.lipids.upper()} and name C2* and not name C2'

    if args.nearest:
        sat += f' and same resid as around {args.distance} resname {args.refmol}'
        unsat += f' and same resid as around {args.distance} resname {args.refmol}'

    return sat, unsat


def process_frames(u, args):
    """
        Calculate the range of frames based on the provided arguments.

        Args:
            args (Namespace): Command-line arguments containing `b`, `e`, `ts`, `bf`, `ef`, and `step`.
            u (mda.Universe): MDAnalysis Universe object

        Returns:
            list: A list of frame indices based on the calculated range.
    """

    time = u.trajectory.totaltime
    total_frames = len(u.trajectory)

    if args.b is not None or args.e is not None or args.ts is not None:
        frames_range = range(int(args.b * total_frames / time) if args.b is not None else 0,
            int(args.e * total_frames / time) if args.e is not None else total_frames - 1,
            int(args.ts * total_frames / time) if args.ts is not None else 1)
    else:
        frames_range = range(args.bf if args.bf is not None else 0,
                             args.ef if args.ef is not None else total_frames,
                             args.step if args.step is not None else 1)

    frames = list(frames_range)

    return frames


def calculate_order_parameters(universe, frames, selection):
    """
        Calculate lipid order parameters for selected frames.

        Args:
            universe: MDAnalysis Universe object
            frames: Frames to analyze
            selection: List of two atom selections for lipid tails

        Returns:
            numpy.ndarray: Order parameters array with shape
            (tails, carbon atoms)
    """
    num_frames = len(universe.trajectory[frames])

    # Pre-compute lipid information
    sat_tail_atoms = universe.select_atoms(selection[0])
    unsat_tail_atoms = universe.select_atoms(selection[1])

    max_sat_num = max(int(atom.name[2:]) for atom in sat_tail_atoms)
    max_unsat_num = max(int(atom.name[2:]) for atom in unsat_tail_atoms)
    num_lipids = len(sat_tail_atoms.split('residue'))
    num_atoms = max(max_unsat_num, max_sat_num)

    # Initialize output array (2 for number of lipid tails, 3 for hydrogens)
    order_params = np.full(
        (num_frames, 2, num_lipids, num_atoms - 1, 3),
        np.nan
    )

    with tqdm(total=num_frames, desc="Processing Frames") as progress_bar:
        for frame_idx, _ in enumerate(universe.trajectory[frames]):
            for tail_idx, tail_selection in enumerate(selection):
                lipids = universe.select_atoms(tail_selection).split('residue')

                for lipid_idx, lipid in enumerate(lipids):
                    for atom in lipid.atoms:
                        atom_num = int(atom.name[2:]) - 2

                        for hydrogen_idx, hydrogen in enumerate(
                                h for h in atom.bonded_atoms if 'H' in h.name
                        ):
                            # Calculate C-H vector and order parameter
                            ch_vector = hydrogen.position - atom.position
                            norm_vector = np.linalg.norm(ch_vector)
                            cos_theta = ch_vector[2] / norm_vector  # [0,0,1] dot product simplifies to z-component
                            order_param = 0.5 * (3 * cos_theta ** 2 - 1)

                            order_params[frame_idx, tail_idx, lipid_idx, atom_num, hydrogen_idx] = order_param

            progress_bar.update(1)

    return order_params


def compute_order_param_stats(order_params):
    """Compute mean and SEM of order parameters after hierarchical averaging."""
    # Step 1: Average over hydrogens (axis=4)
    mean_h = np.nanmean(order_params, axis=4)
    n_h = np.sum(~np.isnan(order_params), axis=4)
    var_h = np.nanvar(order_params, axis=4)
    sem_h = np.sqrt(var_h / n_h)

    # Step 2: Average over lipids (axis=2)
    mean_lipids = np.nanmean(mean_h, axis=2)
    n_lipids = np.sum(~np.isnan(mean_h), axis=2)
    sem_lipids = np.sqrt(np.nansum(sem_h ** 2, axis=2) / (n_lipids ** 2))

    # Step 3: Average over frames (axis=0)
    mean_frames = np.nanmean(mean_lipids, axis=0)
    n_frames = np.sum(~np.isnan(mean_lipids), axis=0)
    sem_frames = np.sqrt(np.nansum(sem_lipids ** 2, axis=0) / (n_frames ** 2))

    # Step 4: Average over frames (axis=0)
    if mean_frames[0, ~np.isnan(mean_frames[0])].shape == mean_frames[1, ~np.isnan(mean_frames[1])].shape:
        mean_tails = np.nanmean(mean_frames, axis=0)
        n_tails = np.sum(~np.isnan(mean_frames), axis=0)
        sem_tails = np.sqrt(np.nansum(sem_frames ** 2, axis=0) / (n_tails ** 2))
        return np.abs(mean_tails), sem_tails

    return np.abs(mean_frames), sem_frames


def block_average(order_params, n_blocks):
    """Compute mean and SEM using block averaging."""
    data = np.nanmean(order_params, axis=(2, 4))
    n_frames = data.shape[0]
    block_size = n_frames // n_blocks
    blocks = []

    for i in range(n_blocks):
        start = i * block_size
        end = (i + 1) * block_size
        block_mean = np.nanmean(data[start:end], axis=0)
        blocks.append(block_mean)

    blocks = np.array(blocks)  # Shape: (n_blocks, ...)
    block_means = np.nanmean(blocks, axis=0)  # Global mean
    block_sem = np.nanstd(blocks, axis=0) / np.sqrt(n_blocks)  # SEM across blocks

    return blocks, block_means, block_sem


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

    parser.add_argument('-odump', type=str,
                        help='Output file for raw binary data, can only be used with option dump (default: %(default)s)',
                        default='order_data.npz',
                        metavar='<.npz>')

    parser.add_argument('-oblocks', type=str,
                        help='Output binary file for blocks averaging data, can only be used with option block-averaging (default: %(default)s)',
                        default='order_blocks.npz',
                        metavar='<.npz>')

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

    parser.add_argument('--dump',
                        help='save unedited data',
                        default=True,
                        action=argparse.BooleanOptionalAction)

    parser.add_argument('--block-averaging', type=int,
                        metavar='<blocks>',
                        default=None,
                        help='Number of blocks for block averaging, if needed (default: 0)')

    args = parser.parse_args()

    u = mda.Universe(args.s, args.f)

    frames = process_frames(u, args)
    selection = select_lipids(args)

    order_params = calculate_order_parameters(u, frames, selection)
    order_params_mean = compute_order_param_stats(order_params)

    if args.dump:
        np.savez(args.odump, order_params=order_params)

    if args.block_averaging is not None:
        blocks, block_means, block_sem = block_average(order_params, args.block_averaging)
        print(blocks)
        np.savez("lipid_order_params.npz", blocks=blocks, means=block_means, sem=block_sem)

    if args.lipids.upper() == 'POPC':
        header = f'Saturated Unsaturated SEM_s SEM_us'
    else:
        header = f'OP SEM'

    np.savetxt(args.o, np.row_stack(order_params_mean).T, header=header, comments='')
