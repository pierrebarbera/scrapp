#!/usr/bin/env python2

import argparse
import os
import sys
import glob

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import scripts.util as util

import subprocess as sub

FNULL = open(os.devnull, 'wb')

# a wrapper script to parallelize many calls to msa_bootstrap, using threading or MPI, depending on how it was called

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
        description="Wrapper script that calls msa_bootstrap in parallel, using threading or MPI, on a given set of directories"
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
        '--num-replicates',
        help="Number of bootstrap replicates to generate.",
        action='store',
        dest='num_reps',
        default=20,
        type=int
    )

    parser.add_argument(
        "--verbose",
        help="Increase output verbosity.",
        action="store_true"
    )

    parser.add_argument(
        "--no-cleanup",
        help="Specify if all intermediate files should be kept (potentially thousands!).",
        action="store_false",
        dest='cleanup'
    )

    parser.add_argument(
        "--model",
        help="Model to use.",
        action="store",
        dest="model",
        type=str
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

def run_func( edge_dir, args ):
    # TODO to minimize concurrent disk footprint: generate single BS MSA just before evaluation, using raxml-ng --bsmsa

    bs_reps_out_dir = os.path.join( args.work_dir, edge_dir, "bs_rep_msas" )
    bs_reps_chk_file = os.path.join( bs_reps_out_dir, "bs_reps_cmd.txt" )
    bs_reps_out_file = os.path.join( bs_reps_out_dir, "bs_reps_log.txt" )

    paths = util.subprogram_commands()

    msa = os.path.join( args.work_dir, edge_dir, "aln.fasta" )

    bs_reps_cmd = [
        paths[ "msa_bootstrap" ],
        msa,
        bs_reps_out_dir,
        str( args.num_reps )
    ]

    if ( not util.call_with_check_file(
        bs_reps_cmd,
        bs_reps_chk_file,
        out_file_path=bs_reps_out_file,
        err_file_path=bs_reps_out_file,
        verbose=args.verbose
    ) ):
        raise RuntimeError( "get_all_bs_reps has failed! (" + bs_reps_out_file + ")" )

    trees_eval_out_dir = os.path.join( args.work_dir, edge_dir, "trees")
    trees_eval_out_file = os.path.join( trees_eval_out_dir, "bs_reps_log.txt" )

    best_trees = glob.glob( os.path.join( args.work_dir, edge_dir, "search", "*.raxml.bestTree" ) )

    if not best_trees:
        raise Exception("No trees were found in `" + os.path.join( args.work_dir, edge_dir, "search") + "` during run_bootstraps_processes" )

    bsrep_msas = glob.glob( os.path.join( bs_reps_out_dir, "replicate_*.fasta" ) )

    for rep_msa in bsrep_msas:
        filename = rep_msa.split( "/" )[-1]
        name = filename.split( "." )[0]

        trees_eval_chk_file = os.path.join( trees_eval_out_dir, name + "_eval_cmd.txt" )

        # reevaluate the tree under this bootstrapped MSA
        trees_eval_cmd = [
            paths[ "raxml-ng" ],
            "--evaluate",
            "--msa", rep_msa,
            "--tree", best_trees[0],
            "--prefix", trees_eval_out_dir + "/" + name,
            "--model", args.model,
            "--threads", "1"
        ]

        if ( not util.call_with_check_file(
            trees_eval_cmd,
            trees_eval_chk_file,
            out_file_path=trees_eval_out_file,
            err_file_path=trees_eval_out_file,
            verbose=args.verbose
        ) ):
            raise RuntimeError( "trees_eval has failed!" )

        # clean up BS MSA
        if args.cleanup:
            os.remove( rep_msa )

        # rename the resulting tree
        sub.call(["ln", "-s", os.path.join( trees_eval_out_dir, name + ".raxml.bestTree" ),
                                os.path.abspath(os.path.join( trees_eval_out_dir, name + ".newick"))], stdout=FNULL)

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
                raise RuntimeError( "msa_bootstrap has failed! "  + f.result() )
    except ImportError:
        from concurrent.futures import ProcessPoolExecutor
        from multiprocessing import Pool
        executor = ProcessPoolExecutor(max_workers=threads)
        pool = Pool(processes=threads)

        results = [pool.apply_async( run_func, args=(edge_dir, args)) for edge_dir in args.edge_dirs]
        for result in results:
            if result.get() > 0:
                raise RuntimeError( "msa_bootstrap has failed!" )

