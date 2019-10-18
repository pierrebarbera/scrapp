#!/usr/bin/env python2

import argparse
import os
import sys
import glob

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import scripts.util as util

# a wrapper script to parallelize many calls to mptp, using threading or MPI, depending on how it was called

def command_line_args_parser():
    """
    Return an instance of argparse that can be used to process command line arguemnts.
    """

    def unit_interval(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
        return x

    # Init an args parser, with a group of required named arguments. It is just nicer to use named
    # arguments than having to rely on their order (i.e., use positional arguments instead).
    parser = argparse.ArgumentParser(
        description="Wrapper script that calls mPTP in parallel, using threading or MPI, on a given set of directories"
    )
    parser_required = parser.add_argument_group('required arguments')

    parser_required.add_argument(
        '--work-dir',
        help="Working directory.",
        action='store',
        dest='work_dir',
        type=str,
        required=True
    )

    parser_required.add_argument(
        'edge_dirs',
        help="Directories to process.",
        action='store',
        type=str,
        nargs='+'
    )

    parser.add_argument(
        '--seed',
        help="Seed for the RNG.",
        action='store',
        dest='seed',
        type=int
    )

    parser.add_argument(
        '--threads',
        help="Number of threads to use if we are not called via mpiexec (there the mpiexec -n determines number of ranks).",
        action='store',
        dest='threads',
        default=1,
        type=int
    )

    parser.add_argument(
        "--verbose",
        help="Increase output verbosity.",
        action="store_true"
    )

    return parser

def command_line_args_postprocessor( args ):
    # Make sure that all paths are fully resolved and dirs have no trailing slashes.
    # args.jplace_file = os.path.abspath( os.path.realpath( args.jplace_file ))

    return args

def command_line_args():
    """
    Return a parsed and processed list of the command line arguments that were provided when
    running this script.
    """

    # Parse the given arguments from the command line, post-process them, return the result.
    parser = command_line_args_parser()
    args = parser.parse_args()
    args = command_line_args_postprocessor( args )
    return args

def is_master():
    try:
        from mpi4py import MPI
        return MPI.COMM_WORLD.Get_rank() == 0
    except ImportError:
        return True

def run_func( edge_dir, args ):
    assert os.path.exists( edge_dir )
    assert os.path.exists( args.work_dir )

    trees = glob.glob( os.path.join( args.work_dir, edge_dir, "trees", "*.newick") )
    for tree in trees:
        # get the name of the file to use as a subdirectory name
        tree_name = os.path.basename( tree ).split(".", 2)[0]

        # set the paths/files accordingly
        mptp_out_dir = os.path.join( args.work_dir, edge_dir, "delimit", tree_name)
        mptp_chk_file = os.path.join( mptp_out_dir, "mptp_cmd.txt" )
        mptp_out_file = os.path.join( mptp_out_dir, "mptp_log.txt" )

        paths = util.subprogram_commands()

        mptp_cmd = [
            paths[ "mptp" ],
            "--tree_file", tree,
            "--ml", "--multi",
            "--output_file",  os.path.join( mptp_out_dir, "mptp_result" )
        ]

        if args.seed:
            mptp_cmd.extend(["--seed", str(args.seed)])

        # do a delimitation per possible rooting

        if ( not util.call_with_check_file(
            mptp_cmd,
            mptp_chk_file,
            out_file_path=mptp_out_file,
            err_file_path=mptp_out_file,
            verbose=args.verbose
        ) ):
            raise RuntimeError( "mptp has failed! (log: " + mptp_out_file + ")" )

    return 0

if __name__ == "__main__":
    args = command_line_args()

    threads = args.threads
    if threads < 1:
        raise RuntimeError( "Invalid number of threads: {}".format(threads) )

    # set up the execution pool, use MPI if available and called with mpiexec
    try:
        from mpi4py import MPI
        from mpi4py.futures import MPIPoolExecutor
        print "setting up MPIPoolExecutor, size ",threads
        executor = MPIPoolExecutor(max_workers=threads)
        futures = []
        for edge_dir in args.edge_dirs:
            futures.append( executor.submit( run_func, edge_dir, args ) )

        # wait for all processes to return
        for f in futures:
            if f.result() > 0:
                raise RuntimeError( "mptp has failed! " + f.result() )
    except ImportError:
        from concurrent.futures import ProcessPoolExecutor
        from multiprocessing import Pool
        executor = ProcessPoolExecutor(max_workers=threads)
        pool = Pool(processes=threads)
        results = [pool.apply_async( run_func, args=(edge_dir, args)) for edge_dir in args.edge_dirs]
        for result in results:
            if result.get() > 0:
                raise RuntimeError( "mptp has failed!" )

