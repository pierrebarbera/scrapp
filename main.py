#!/usr/bin/python

import argparse
import sys
import os
import subprocess

# from parallel_decorators.parallel_decorators import vectorize_parallel, is_master

# ==================================================================================================
#     Sub-Program Commands
# ==================================================================================================

def subprogram_commands():
    """
    Return a list of the sub program commands.
    """

    # Get the dir where the script is located. The first variant expects that this is the main
    # script. The second one might however break on a different OS. Needs to be tested.
    # basedir = os.path.abspath( os.path.dirname(sys.argv[0]) )
    basedir = os.path.dirname( os.path.abspath( os.path.realpath( __file__ )))

    # We currently expect all sub-programs to be located in sub-directories of our tool.
    # Maybe this should go into a config file which also allows to use actual commands,
    # in case that raxml-ng or mptp is actually installed on the system already.
    paths = {
        "alignment_splitter" : basedir + "/genesis/bin/apps/",
        "mptp"               : basedir + "/mptp/bin/mptp",
        "raxml-ng"           : basedir + "/raxml-ng/..."
    }
    return paths

def subprograms_exists( paths ):
    """
    Check whether all subprograms actually exists. Throw otherwise.
    """

    for name, cmd in paths.iteritems():
        if not os.path.isfile( cmd ):
            raise RuntimeError(
                "Subprogram '" + name + "' not found at '" + cmd + "'. Please run setup first."
            )

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
        action='store', dest='jplace_file',
        required=True
    )
    parser_required_named_arg_group.add_argument(
        "-a", "--alignment",
        help="The alignment file path containing the alignment of the reads. Fasta or Phylip format.",
        action='store', dest='aln_file',
        required=True
    )

    # Add optional args.
    parser.add_argument(
        '-o', '--output',
        help="Output file path for the Newick file containg species counts per branch.",
        action='store', dest='output_file',
        default=None
    )
    parser.add_argument(
        '-w', '--work-dir',
        help="Directory path for intermediate work files.",
        action='store', dest='work_dir',
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
        default=1.0
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

def call_with_check_file( cmd_to_call, check_file_path, verbose=False ):
    """
    Run a shell command if it was not run before, using a check file to find out whether we ran it
    before or not.
    """

    if verbose:
        print "Running command: " + cmd_to_call

    # If the check file exists, it has to contain the exact same command.
    if os.path.isfile( check_file_path ):
        with open( check_file_path, 'r') as check_file_handle:
            check_file_content = check_file_handle.read()

        if check_file_content == cmd_to_call:
            print "Already did this step. Skipping."
        else:
            raise RuntimeError(
                "Check file '" + check_file_path + "' already exists but has unexpected content. "
                "This most likely means that you ran SCRAPP before, using the same work "
                "directory, but different input files."
            )

    # If the checkfile does not exist, run the command, then create the checkfile and write
    # the command to it.
    else:
        # subprocess.call( cmd_to_call, shell=True );

        if not os.path.exists( os.path.dirname( check_file_path )):
            os.makedirs( os.path.dirname( check_file_path ))
        with open( check_file_path, "w") as check_file_handle:
            check_file_handle.write( cmd_to_call )

# ==================================================================================================
#     Main Function
# ==================================================================================================

if __name__ == "__main__":
    paths = subprogram_commands()
    args  = command_line_args()

    # Print some verbose output about args and params etc.
    if args.verbose:
        print "Command line arguments:", str(args)[len("Namespace("):-1]
        print "Subprogram paths:", paths

    # Check whether all sub programs exist.
    subprograms_exists( paths )

    # Create the work dir to store our stuff.
    if not os.path.exists( args.work_dir ):
        os.makedirs( args.work_dir )

    # Call Genesis to split Jplace file into alignments per branch.
    print "Splitting alignment into per-branch alignments using jplace placements."
    alignment_splitter_check_file = os.path.join( args.work_dir, "alignment_splitter_check.txt" )
    alignment_splitter_cmd = " ".join([
        paths[ "alignment_splitter" ],
        args.jplace_file,
        args.aln_file
    ])
    call_with_check_file( alignment_splitter_cmd, alignment_splitter_check_file, args.verbose )
