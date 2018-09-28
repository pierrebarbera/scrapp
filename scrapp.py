#!/usr/bin/python

import argparse
import glob
import os
import subprocess
import sys
import scripts.util as util
import time
import pprint

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

    parser.add_argument(
        '--min-queries',
        help="If an edge contains a number of unique queries below this value, ignore the edge.",
        action='store',
        dest='min_query',
        type=int,
        default=4
    )
    parser.add_argument(
        '-c', '--cluster-above',
        help="If an edge contains a number of unique queries above this value, apply OTU clustering.",
        action='store',
        dest='max_query',
        type=int,
        default=500
    )
    # parser.add_argument(
    #     '-p', '--parallelization',
    #     help="Parallelization strategy to use. Either 'threads' or 'mpi'.",
    #     action='store', dest='parallelization',
    #     choices=[ "threads", "mpi" ],
    #     default="threads"
    # )
    # parser.add_argument(
    #     '-o', '--out-dir',
    #     help="Output directory.",
    #     action='store',
    #     dest='out_dir',
    #     type=str,
    #     default=base_dir_
    # )
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
    # if args.output_file is None:
    #     # Get file name and extension and append/replace depending on the extension.
    #     jfn, jfe = os.path.splitext( args.jplace_file )
    #     if jfe == ".jplace":
    #         args.output_file = jfn + ".newick"
    #     else:
    #         args.output_file = args.jplace_file + ".newick"

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
    # args.output_file = os.path.abspath( os.path.realpath( args.output_file ))
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

def mkdirp( path ):
    if not os.path.exists( path ):
        os.mkdir( path )

# ==================================================================================================
#     Main Function
# ==================================================================================================

if __name__ == "__main__":
    pp = pprint.PrettyPrinter(indent=4)
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
    # We only do that once in the master rank.
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
            str(args.min_weight),
            str(args.min_query),
            str(args.max_query)
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
    #     OTU Clustering of queries
    # -------------------------------------------------------------------------
    @vectorize_parallel( method = 'adaptive', num_procs = 1 )
    def run_swarm_processes( edge_dir, work_dir ):
        swarm_out_dir = os.path.join( args.work_dir, edge_dir, "swarm")
        swarm_chk_file = os.path.join( swarm_out_dir, "swarm_cmd.txt" )
        swarm_out_file = os.path.join( swarm_out_dir, "swarm_log.txt" )

        sequences = os.path.join( args.work_dir, edge_dir, "aln.fasta" )
        stripped_sequences = os.path.join( args.work_dir, edge_dir, "stripped.fasta" )

        # only proceed if stripped.fasta exists
        if (not os.path.isfile(stripped_sequences)):
            return 0

        otu_path = os.path.join( swarm_out_dir, "otus.fasta")

        swarm_cmd = [
            paths[ "swarm" ],
            "--threads", "1",
            "--fastidious",
            "-w", otu_path,
            stripped_sequences
        ]

        if ( not call_with_check_file(
            swarm_cmd,
            swarm_chk_file,
            out_file_path=swarm_out_file,
            err_file_path=swarm_out_file,
            verbose=args.verbose
        ) ):
            raise RuntimeError( "swarm has failed!" )

        # --------------------------------------------------
        #   transform the result otus back to aligned ones
        # --------------------------------------------------
        map_back_out_dir = swarm_out_dir
        map_back_chk_file = os.path.join( map_back_out_dir, "map_back_cmd.txt" )
        map_back_out_file = os.path.join( map_back_out_dir, "map_back_log.txt" )

        map_back_cmd = [
            paths[ "otu_map_back" ],
            otu_path,
            sequences,
            os.path.join( args.work_dir, edge_dir, "aligned_otus.fasta")
        ]

        if ( not call_with_check_file(
            map_back_cmd,
            map_back_chk_file,
            out_file_path=map_back_out_file,
            err_file_path=map_back_out_file,
            verbose=args.verbose
        ) ):
            raise RuntimeError( "map_back has failed!" )

        print "done", mpi_rank(), 0
        return 0

    runtime = time.time()
    run_swarm_processes( edge_list, args.work_dir )
    runtime = time.time() - runtime
    runtimes.append({"name":"swarm", "time":str(runtime)})

    # exit()

    # -------------------------------------------------------------------------
    #     RAxML Tree Inferrence
    # -------------------------------------------------------------------------

    # Create a parallel function that either runs on multiple MPI nodes,
    # each of them running one RAxML instance with as many threads as specified in the CLI,
    # or, if we are not using MPI, run the parallel loop single threaded,
    # but use the threads again internally for the RAxML instance.


    if is_master():

        # prepare for pargenes by copying, renaming msa files, into a temp dir
        import shutil
        tmp_dir = os.path.join( args.work_dir, "tmp" )
        mkdirp( tmp_dir )
        for edge_dir in edge_list:
            msa = os.path.join( args.work_dir, edge_dir, "aligned_otus.fasta" )
            if (not os.path.isfile( msa )):
                if (os.path.isfile( os.path.join(args.work_dir, edge_dir, "stripped.fasta") )):
                    raise RuntimeError( "stripped.fasta exists, but there is no otu file, something must have gone wrong!" )
                msa = os.path.join( args.work_dir, edge_dir, "aln.fasta" )
            edge_string = edge_dir.split("/")[-2]
            shutil.copyfile(msa, os.path.join(tmp_dir, edge_string + ".fasta"))

        # call pargenes
        pargenes_chk_file = os.path.join( args.work_dir, "pargenes_cmd.txt" )
        pargenes_out_file = os.path.join( args.work_dir, "pargenes_log.txt" )

        pargenes = os.path.join(base_dir_, "deps/ParGenes/pargenes/pargenes.py")

        tmp_out_dir = os.path.join( args.work_dir, "tmp_out" )
        # os.mkdir(tmp_out_dir)

        pargenes_cmd = ["python", pargenes,
            "--alignments-dir", tmp_dir,
            "--output-dir", tmp_out_dir,
            "--datatype", "nt",
            "--cores", str(args.num_threads),
            "--scheduler", "openmp",
            "--continue",
            "--raxml-global-parameters-string", "--model GTR+G"
        ]

        runtime = time.time()
        if ( not call_with_check_file(
            pargenes_cmd,
            pargenes_chk_file,
            out_file_path=pargenes_out_file,
            err_file_path=pargenes_out_file,
            verbose=args.verbose
        ) ):
            raise RuntimeError( "pargenes has failed!" )
        runtime = time.time() - runtime
        runtimes.append({"name":"pargenes", "time":str(runtime)})

        # copy the results back to their appropriate directories
        for edge_dir in edge_list:
            pargenes_out_dir = os.path.join( args.work_dir, edge_dir, "search/" )
            edge_string = edge_dir.split("/")[-2]
            res = os.path.join( tmp_out_dir, "mlsearch_run/results", edge_string + "_fasta" )
            for filename in glob.glob(os.path.join(res, "*.*")):
                # print "copy ", filename, " to ", pargenes_out_dir
                mkdirp( pargenes_out_dir )
                shutil.copy( filename, pargenes_out_dir )

    # -------------------------------------------------------------------------
    #     get all possible rootings per tree
    # -------------------------------------------------------------------------
    @vectorize_parallel( method = 'adaptive', num_procs = 1 )
    def run_rootings_processes( edge_dir, work_dir ):
        rootings_out_dir = os.path.join( args.work_dir, edge_dir, "trees")
        rootings_chk_file = os.path.join( rootings_out_dir, "rootings_cmd.txt" )
        rootings_out_file = os.path.join( rootings_out_dir, "rootings_log.txt" )

        bestTree = glob.glob( os.path.join( args.work_dir, edge_dir, "search", "*.raxml.bestTree" ) )

        rootings_cmd = [
            paths[ "get_all_rootings" ],
            bestTree[0],
            rootings_out_dir
        ]

        if ( not call_with_check_file(
            rootings_cmd,
            rootings_chk_file,
            out_file_path=rootings_out_file,
            err_file_path=rootings_out_file,
            verbose=args.verbose
        ) ):
            raise RuntimeError( "get_all_rootings has failed!" )

        print "done", mpi_rank(), 0
        return 0

    runtime = time.time()
    run_rootings_processes( edge_list, args.work_dir )
    runtime = time.time() - runtime
    runtimes.append({"name":"get_all_rootings", "time":str(runtime)})

    # -------------------------------------------------------------------------
    #     Species Delimitation
    # -------------------------------------------------------------------------

    @vectorize_parallel( method = 'adaptive', num_procs = 1 )
    def run_mptp_processes( edge_dir, work_dir ):

        # switch this via an option (best tree only vs. all rootings vs. longest edge rooting)
        trees = glob.glob( os.path.join( args.work_dir, edge_dir, "trees", "*.newick") )
        for tree in trees:
            # get the name of the file to use as a subdirectory name
            tree_name = os.path.basename( tree ).split(".", 2)[0]

            # set the paths/files accordingly
            mptp_out_dir = os.path.join( args.work_dir, edge_dir, "delimit", tree_name)
            mptp_chk_file = os.path.join( mptp_out_dir, "mptp_cmd.txt" )
            mptp_out_file = os.path.join( mptp_out_dir, "mptp_log.txt" )

            mptp_cmd = [
                paths[ "mptp" ],
                "--tree_file", tree,
                "--ml", "--multi",
                "--output_file",  os.path.join( mptp_out_dir, "mptp_result" )
            ]

            # do a delimitation per possible rooting

            if ( not call_with_check_file(
                mptp_cmd,
                mptp_chk_file,
                out_file_path=mptp_out_file,
                err_file_path=mptp_out_file,
                verbose=False
            ) ):
                raise RuntimeError( "mptp has failed!" )

        print "done", mpi_rank(), 0
        return 0

    runtime = time.time()
    run_mptp_processes( edge_list, args.work_dir )
    runtime = time.time() - runtime
    runtimes.append({"name":"mptp", "time":str(runtime)})

    # -------------------------------------------------------------------------
    #     Summarize Delimitation Results
    # -------------------------------------------------------------------------
    if is_master():
        import scripts.mptp as mptp
        import scripts.tea as tea

        output = tea.TEA()

        # for all reference edges (that have results)s
        for d in edge_list:
            d = os.path.join(args.work_dir, d, "delimit")

            # for all possible runs/rootings of the delimitation
            file_paths = glob.glob( os.path.join(d, "edge_*/mptp_result.txt" ) )
            res = []
            for path in file_paths:
                # parse the results
                res.append( mptp.parse(path) )
            # summarize them
            summary = mptp.summarize( res )

            # add the summary to the overall result structure
            ref_edge_id = d.split("/")[-3].split("_")[-1]
            output.add_annotation("species-count", ref_edge_id, summary)

        if args.verbose:
            output.to_stream(sys.stdout)
        output.to_file( os.path.join( args.work_dir, "summary.tea" ) )

    if is_master():
        print "Finished!"
        if args.verbose:
            pp.pprint( runtimes )
