#!/usr/bin/python

import argparse
import glob
import os
import subprocess
import sys
import scripts.util as util
import time

base_dir_ = os.path.dirname( os.path.realpath(__file__) )

execfile( os.path.join( base_dir_, "deps/parallel_decorators/parallel_decorators.py" ) )


# ==================================================================================================
#     Command Line Args
# ==================================================================================================

def command_line_args_parser():
    """
    Return an instance of argparse that can be used to process command line arguemnts.
    """

    # Init an args parser, with a group of required named arguments. It is just nicer to use named
    # arguments than having to rely on their order (i.e., use positional arguments instead).
    parser = argparse.ArgumentParser(
        description="Pipeline wrapper script that calculates species counts for each branch of a "
        "reference tree from phylogenetic placement of reads on that tree."
    )
    parser_required_named_arg_group = parser.add_argument_group('required named arguments')

    # Add required named args.
    parser_required_named_arg_group.add_argument(
        "-j", "--jplace",
        help="The jplace file path containing the tree and the placement of reads on its branches.",
        action='store',
        dest='jplace_file',
        type=str,
        required=True
    )
    parser_required_named_arg_group.add_argument(
        "-a", "--alignment",
        help="The alignment file path containing the alignment of the reads. Fasta or Phylip format.",
        action='store',
        dest='aln_file',
        type=str,
        required=True
    )

    # Add optional args.
    parser.add_argument(
        '-t', '--num-threads',
        help="Number of threads to run for parallelization.",
        action='store',
        dest='num_threads',
        type=int,
        default=0
    )
    # parser.add_argument(
    #     '-p', '--parallelization',
    #     help="Parallelization strategy to use. Either 'threads' or 'mpi'.",
    #     action='store', dest='parallelization',
    #     choices=[ "threads", "mpi" ],
    #     default="threads"
    # )
    parser.add_argument(
        '-o', '--output',
        help="Output file path for the Newick file containg species counts per branch.",
        action='store',
        dest='output_file',
        type=str,
        default=None
    )
    parser.add_argument(
        '-w', '--work-dir',
        help="Directory path for intermediate work files.",
        action='store',
        dest='work_dir',
        type=str,
        default="work"
    )
    parser.add_argument(
        "--verbose",
        help="Increase output verbosity.",
        action="store_true"
    )

    # Add min weight arg, restricted to a certain range, also optional.
    def min_weight_float(x):
        x = float(x)
        if x <= 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range (0.0, 1.0]"%(x,))
        return x
    parser.add_argument(
        '--min-weight',
        help="Minimum weight threshold for placements. Everything below is filtered out.",
        action='store', dest='min_weight',
        type=min_weight_float,
        default=0.5
    )

    return parser

def command_line_args_postprocessor( args ):
    """
    Given the result of argsparse.parse_args(), this function does some specific post-processing
    that we want for our command line arguments.
    """

    # If the user did not specify an output file, use the name of the Jplace file, appending
    # or replacing the the file extension to Newick.
    if args.output_file is None:
        # Get file name and extension and append/replace depending on the extension.
        jfn, jfe = os.path.splitext( args.jplace_file )
        if jfe == ".jplace":
            args.output_file = jfn + ".newick"
        else:
            args.output_file = args.jplace_file + ".newick"

    # If user did not provide number of threads, use all available ones.
    if args.num_threads == 0:
        import multiprocessing
        args.num_threads = multiprocessing.cpu_count()

    # Translate parallelization method to the name used by the decorators.
    # if args.parallelization.lower() == "threads":
    #     args.parallelization = "processes"
    # elif args.parallelization.lower() == "mpi":
    #     args.parallelization = "MPI"
    # else:
    #     raise RuntimeError( "Invalid parallelization method: '" + args.parallelization + "'." )

    # Make sure that all paths are fully resolved and dirs have no trailing slashes.
    args.jplace_file = os.path.abspath( os.path.realpath( args.jplace_file ))
    args.aln_file    = os.path.abspath( os.path.realpath( args.aln_file ))
    args.output_file = os.path.abspath( os.path.realpath( args.output_file ))
    args.work_dir    = os.path.abspath( os.path.realpath( args.work_dir ))

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

# ==================================================================================================
#     Helper Functions
# ==================================================================================================

def call_with_check_file(
    cmd_to_call,
    check_file_path,
    out_file_path=None,
    err_file_path=None,
    verbose=False
):
    """
    Run a shell command if it was not run before, using a check file to find out whether we ran it
    before or not.
    """

    cmd_string = " ".join(cmd_to_call)

    if verbose:
        print "Running command: " + cmd_string

    # If the check file exists, it has to contain the exact same command.
    if os.path.isfile( check_file_path ):
        with open( check_file_path, 'r') as check_file_handle:
            check_file_content = check_file_handle.read()

        if check_file_content == cmd_string:
            if verbose:
                print "Already did this step. Skipping."
            return True
        else:
            raise RuntimeError(
                "Check file '" + check_file_path + "' already exists but has unexpected content. "
                "This most likely means that you ran SCRAPP before, using the same work "
                "directory, but different input files."
            )

    # If the checkfile does not exist, run the command, then create the checkfile and write
    # the command to it.
    else:
        # Prepare out and err files, if needed.
        out_file = None
        if out_file_path is not None:
            if not os.path.exists( os.path.dirname( out_file_path )):
                os.makedirs( os.path.dirname( out_file_path ))
            out_file = open( out_file_path, "w" )
        err_file = None
        if err_file_path is not None:
            if not os.path.exists( os.path.dirname( err_file_path )):
                os.makedirs( os.path.dirname( err_file_path ))
            err_file = open( err_file_path, "w" )

        # Call the command and record its exit code.
        success = ( subprocess.call( cmd_to_call, stdout=out_file, stderr=err_file ) == 0 )

        # If we were not successfull, end the function here.
        if not success:
            return success

        # Only if the command returned successfully, create the checkfile.
        if not os.path.exists( os.path.dirname( check_file_path )):
            os.makedirs( os.path.dirname( check_file_path ))
        with open( check_file_path, "w") as check_file_handle:
            check_file_handle.write( cmd_string )
        return success

# ==================================================================================================
#     Main Function
# ==================================================================================================

if __name__ == "__main__":
    # Get all needed input.
    paths = util.subprogram_commands()
    args  = command_line_args()

    runtimes = []

    # -------------------------------------------------------------------------
    #     Initial Master Work
    # -------------------------------------------------------------------------

    if is_master():
        print "Running SCRAPP"

        # Print some verbose output about args and params etc.
        if args.verbose:
            print "Command line arguments:", str(args)[len("Namespace("):-1]
            print "Subprogram paths:", paths

        # Check whether all sub programs exist.
        # subprograms_exist( paths )

        # Create the work dir to store our stuff.
        if not os.path.exists( args.work_dir ):
            os.makedirs( args.work_dir )

    # -------------------------------------------------------------------------
    #     Split Alignment
    # -------------------------------------------------------------------------

    # Call Genesis to split Jplace file into alignments per branch.
    # We only do that onces in the master rank.
    if is_master():
        print "Splitting alignment into per-branch alignments using jplace placements."

        # Compose the command line args for the call, then execute it.
        aln_splitter_chk_file = os.path.join( args.work_dir, "alignment_splitter_cmd.txt" )
        aln_splitter_out_file = os.path.join( args.work_dir, "alignment_splitter_log.txt" )

        aln_splitter_cmd = [
            paths[ "alignment_splitter" ],
            args.jplace_file,
            args.aln_file,
            args.work_dir,
            str(args.min_weight)
        ]
        runtime = time.time()
        succ = call_with_check_file(
            aln_splitter_cmd,
            aln_splitter_chk_file,
            out_file_path=aln_splitter_out_file,
            err_file_path=aln_splitter_out_file,
            verbose=args.verbose
        )
        runtime = time.time() - runtime
        runtimes.append({"name":"alignment_splitter", "time":str(runtime)})

        # We only continue with the script if the alignment splitting was successfull.
        if not succ:
            print "Could not split the alignment. See log file for details:", aln_splitter_out_file
            sys.exit(1)

        # The result of alignment splitting is stored in sub-directories in our work dir.
        # The list of those dirs is what we need to process now.
        edge_list = glob.glob( args.work_dir + "/edge_*/" )

        # User output
        print "Processing", len(edge_list), "edges."
        if args.verbose:
            for edge in edge_list:
                print "  - " + edge

    else:
        # For non-master ranks, we create a dummy list, which is passed to the parallel function.
        # This is then internally overriden by a broadcast of the actual list of the master rank.
        edge_list = []

    # edge_list = [edge_list[0]]

    # -------------------------------------------------------------------------
    #     RAxML Tree Inferrence
    # -------------------------------------------------------------------------

    # Create a parallel function that either runs on multiple MPI nodes,
    # each of them running one RAxML instance with as many threads as specified in the CLI,
    # or, if we are not using MPI, run the parallel loop single threaded,
    # but use the threads again internally for the RAxML instance.
    @vectorize_parallel( method = 'adaptive', num_procs = 1 )
    def run_raxml_processes( edge_dir, work_dir ):
        raxml_out_dir = os.path.join( args.work_dir, edge_dir, "search")
        raxml_chk_file = os.path.join( raxml_out_dir, "raxml_cmd.txt" )
        raxml_out_file = os.path.join( raxml_out_dir, "raxml_log.txt" )

        raxml_cmd = [
            paths[ "raxml-ng" ], "--search",
            "--msa", os.path.join( args.work_dir, edge_dir, "aln.phylip" ),
            "--threads", str(args.num_threads),
            "--prefix", raxml_out_dir + "/s",
            "--model", "GTR+G"
        ]

        if ( not call_with_check_file(
            raxml_cmd,
            raxml_chk_file,
            out_file_path=raxml_out_file,
            err_file_path=raxml_out_file,
            verbose=args.verbose
        ) ):
            raise RuntimeError( "raxml has failed!" )


        print "done", mpi_rank(), succ
        return succ

    runtime = time.time()
    run_raxml_processes( edge_list, args.work_dir )
    runtime = time.time() - runtime
    runtimes.append({"name":"raxml-ng", "time":str(runtime)})

    # Tests
    # @vectorize_parallel( method = 'MPI' )
    # def simple_test( list ):
    #     print list, mpi_rank()
    # simple_test( [1,2,3] )
    # print "whoami", mpi_rank(), mpi_size(), is_master()

    # -------------------------------------------------------------------------
    #     Species Delimitation
    # -------------------------------------------------------------------------
    @vectorize_parallel( method = 'adaptive', num_procs = 1 )
    def run_mptp_processes( edge_dir, work_dir ):
        mptp_out_dir = os.path.join( args.work_dir, edge_dir, "delimit")
        mptp_chk_file = os.path.join( mptp_out_dir, "mptp_cmd.txt" )
        mptp_out_file = os.path.join( mptp_out_dir, "mptp_log.txt" )

        mptp_cmd = [
            paths[ "mptp" ],
            "--tree_file", os.path.join( args.work_dir, edge_dir, "search", "s.raxml.bestTree" ), #TODO actual name
            "--ml", "--multi",
            "--output_file",  os.path.join(mptp_out_dir, "mptp_result.txt")
        ]

        if ( not call_with_check_file(
            mptp_cmd,
            mptp_chk_file,
            out_file_path=mptp_out_file,
            err_file_path=mptp_out_file,
            verbose=args.verbose
        ) ):
            raise RuntimeError( "mptp has failed!" )


        print "done", mpi_rank(), succ
        return succ

    runtime = time.time()
    run_mptp_processes( edge_list, args.work_dir )
    runtime = time.time() - runtime
    runtimes.append({"name":"mptp", "time":str(runtime)})

    if is_master():
        print "Finished!"
        print runtimes
