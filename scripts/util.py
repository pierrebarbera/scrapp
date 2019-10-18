import os
import platform
import subprocess as sub
import urllib2
import time
from collections import OrderedDict

scripts_dir_ = os.path.dirname( os.path.realpath(__file__) )

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
            out_file = open( out_file_path, "w+" )
        err_file = None
        if err_file_path is not None:
            if not os.path.exists( os.path.dirname( err_file_path )):
                os.makedirs( os.path.dirname( err_file_path ))
            err_file = open( err_file_path, "w+" )

        # Call the command and record its exit code.
        success = ( sub.call( cmd_to_call, stdout=out_file, stderr=err_file ) == 0 )

        # If we were not successfull, end the function here.
        if not success:
            return success

        # Only if the command returned successfully, create the checkfile.
        if not os.path.exists( os.path.dirname( check_file_path )):
            os.makedirs( os.path.dirname( check_file_path ))
        with open( check_file_path, "w+") as check_file_handle:
            check_file_handle.write( cmd_string )
        return success

def call_wrapped( name, edge_list, args, extra=[] ):
    cmd = ["mpiexec", str(args.num_threads)] if args.parallel == "mpi" else []

    wrapper = os.path.join(scripts_dir_, "{}_wrap.py".format( name ))
    assert os.path.isfile( wrapper )
    assert os.access( wrapper, os.X_OK )

    cmd.extend([ wrapper, "--work-dir", args.work_dir ])
    if args.parallel == "threads":
        cmd.extend(["--threads", str(args.num_threads)])
    if args.seed:
        cmd.extend(["--seed", str(args.seed)])
    if args.verbose:
        cmd.extend(["--verbose"])
    if extra:
        cmd.extend( extra )
    cmd.extend( edge_list )

    chk_file = os.path.join( args.work_dir, "{}_wrap_cmd.txt".format( name ) )
    out_file = os.path.join( args.work_dir, "{}_wrap_log.txt".format( name ) )

    if args.verbose: print "running {}!".format( name )
    runtime = time.time()
    if ( not call_with_check_file(
        cmd,
        chk_file,
        out_file_path=out_file,
        err_file_path=out_file,
        verbose=args.verbose
    ) ):
        raise RuntimeError( "{}_wrap has failed!".format( name ) )
    runtime = time.time() - runtime
    return {"name":name, "time":runtime}

# ==================================================================================================
#     globals
# ==================================================================================================

# Get the dir where the script is located. The first variant expects that this is the main
# script. The second one might however break on a different OS. Needs to be tested.
# basedir = os.path.abspath( os.path.dirname(sys.argv[0]) )
basedir = os.path.abspath( os.path.realpath( os.path.join(
            os.path.dirname( os.path.abspath( os.path.realpath( __file__ ))),
            "../deps/")))

# We currently expect all sub-programs to be located in sub-directories of our tool.
# Maybe this should go into a config file which also allows to use actual commands,
# in case that raxml-ng or mptp is actually installed on the system already.
prog_paths = OrderedDict([
    ("pargenes",           [basedir + "/ParGenes/is_installed_"]),
    ("raxml-ng",           [basedir + "/ParGenes/raxml-ng/bin/"]),
    ("libgenesis.so",      [basedir + "/genesis/bin/"]),
    ("alignment_splitter", [basedir + "/genesis/bin/apps/"]),
    ("get_rooting",        [basedir + "/genesis/bin/apps/"]),
    ("phy2fasta",          [basedir + "/genesis/bin/apps/"]),
    ("msa_bootstrap",      [basedir + "/genesis/bin/apps/"]),
    ("mptp",               [basedir + "/mptp/bin/"])
])

FNULL = open(os.devnull, 'wb')

# ==================================================================================================
#     general util
# ==================================================================================================

def get_platform():
    """
    Return a descriptive tuple specifying platform info.
    """
    oper = platform.system().lower()
    arch = platform.machine()

    intrin = []
    tested_vecs = ["sse3", "avx", "avx2", "avx512"]

    for vec in tested_vecs:
        if sub.call(["grep", vec, "/proc/cpuinfo"], stdout=FNULL, stderr=FNULL) == 0:
            intrin = [vec] + intrin

    return (oper, arch, intrin)

def try_fetch(src_url, target_path):
    """
    Tries to download file at <src_url> to <target_path>.
    Returns file path on success, None otherwise.
    """
    written_filepath = ""
    try:
        f = urllib2.urlopen( src_url )
        print "downloading " + src_url

        if not os.path.exists(target_path):
            os.makedirs(target_path)

        # Open our local file for writing
        with open(os.path.join(target_path, os.path.basename(src_url)), "wb+") as local_file:
            written_filepath =  local_file.name
            print "... to " + written_filepath
            local_file.write(f.read())

    except urllib2.HTTPError, e:
        written_filepath = None
    except urllib2.URLError, e:
        written_filepath = None

    return written_filepath

def extract(filepath):
    """
    Extracts an archive based on its ending.
    """
    name, ending = os.path.splitext( filepath )
    if ending == ".zip":
        import zipfile
        zip_ref = zipfile.ZipFile( filepath, 'r' )
        zip_ref.extractall( os.path.dirname( filepath ) )
        zip_ref.close()
    else:
        raise RuntimeError( "Unrecognized extension: " + ending + " for file " + filepath )


def try_resolve_raxmlng(machine = get_platform()):
    """
    Tries to resolve dependency using a cascade of increasingly desperate measures.
    """
    entry = prog_paths["raxml-ng"]
    binpath = entry[0]

    (oper, arch, intrin) = machine

    raxmldir = os.path.join(basedir, "ParGenes/raxml-ng")

    sub.call(["mkdir", "-p", os.path.join(raxmldir, "build")], stdout=FNULL)
    sub.call(["cmake", ".."], cwd=os.path.join(raxmldir, "build"), stdout=FNULL)
    return sub.call(["make"], cwd=os.path.join(raxmldir, "build"), stdout=FNULL)

    # # if binary doesn't exist, see if binary exists online
    # url = entry[1]
    # core_rax = "_" + oper + "_" + arch
    # base_rax_url = url + core_rax + ".zip"
    # # mpi_rax_url = url + core_rax + "_MPI" + ".zip"

    # print "Trying to download raxml-ng binary..."
    # filepath = try_fetch(base_rax_url, binpath)

    # if filepath:
    #     print "Extracting..."
    #     extract(filepath)
    #     print "Done!"

    # TODO: set access modifier

    # TODO: case: cannot fetch by url -> try source


def try_resolve_mptp(machine = get_platform()):
    mptpdir = os.path.join(basedir, "mptp")
    # modify the config file, disabling gsl
    sub.call(["sed", "-i", '45,46 s/^/#/', "configure.ac"], cwd=mptpdir)
    # autogen, configure, make
    sub.call([os.path.join(mptpdir, "autogen.sh")], cwd=mptpdir, stdout=FNULL)
    sub.call([os.path.join(mptpdir, "configure")], cwd=mptpdir, stdout=FNULL)
    return sub.call(["make", "-C", mptpdir], stdout=FNULL)

def try_resolve_genesis(machine = get_platform()):
    genesisdir = os.path.join(basedir, "genesis")
    # make genesis
    return sub.call(["make", "-C", genesisdir], stdout=FNULL)

def try_resolve_alignment_splitter(machine = get_platform()):
    genesisdir = os.path.join(basedir, "genesis")
    # ensure the symlink exists
    sub.call(["ln", "-sft", os.path.join( genesisdir, "apps" ), os.path.abspath(os.path.join(basedir, "../src/alignment_splitter.cpp"))], stdout=FNULL)

    # make update on genesis
    return sub.call(["make", "update", "-C", genesisdir], stdout=FNULL)

def try_resolve_phy2fasta(machine = get_platform()):
    genesisdir = os.path.join(basedir, "genesis")
    # ensure the symlink exists
    sub.call(["ln", "-sft", os.path.join( genesisdir, "apps" ), os.path.abspath(os.path.join(basedir, "../src/phy2fasta.cpp"))], stdout=FNULL)

    # make update on genesis
    return sub.call(["make", "update", "-C", genesisdir], stdout=FNULL)

def try_resolve_msa_bootstrap(machine = get_platform()):
    genesisdir = os.path.join(basedir, "genesis")
    # ensure the symlink exists
    sub.call(["ln", "-sft", os.path.join( genesisdir, "apps" ), os.path.abspath(os.path.join(basedir, "../src/msa_bootstrap.cpp"))], stdout=FNULL)

    # make update on genesis
    return sub.call(["make", "-C", genesisdir], stdout=FNULL)

def try_resolve_pargenes(machine = get_platform()):
  pargenesdir = os.path.join(basedir, "ParGenes")
  # print(pargenesdir)
  res = sub.call(['./install.sh'], cwd=pargenesdir, shell=True, stdout=FNULL)
  if (0 == res):
    sub.call(['touch', os.path.join(pargenesdir, "is_installed_pargenes")])
  return res

def try_resolve_get_rooting(machine = get_platform()):
    genesisdir = os.path.join(basedir, "genesis")
    # ensure the symlink exists
    sub.call(["ln", "-sft", os.path.join( genesisdir, "apps" ), os.path.abspath(os.path.join(basedir, "../src/get_rooting.cpp"))], stdout=FNULL)

    # make update on genesis
    return sub.call(["make", "update", "-C", genesisdir], stdout=FNULL)

def try_resolve(name, machine = get_platform()):
    print "Trying to resolve " + name + "..."
    if name == "raxml-ng":
        return try_resolve_raxmlng( machine )
    elif name == "mptp":
        return try_resolve_mptp( machine )
    elif name == "libgenesis.so":
        return try_resolve_genesis( machine )
    elif name == "alignment_splitter":
        return try_resolve_alignment_splitter( machine )
    elif name == "get_rooting":
        return try_resolve_get_rooting( machine )
    elif name == "phy2fasta":
        return try_resolve_phy2fasta( machine )
    elif name == "msa_bootstrap":
        return try_resolve_msa_bootstrap( machine )
    elif name == "pargenes":
        return try_resolve_pargenes( machine )
    else:
        raise RuntimeError( "No such subprogram name: " + name )


# ==================================================================================================
#     Sub-Program Commands
# ==================================================================================================

def subprogram_commands():
    """
    Return a list of the sub program commands.
    """

    paths = OrderedDict()

    for name, etc in prog_paths.iteritems():
        paths[name] = etc[0] + name

    return paths

def subprograms_exist( paths ):
    """
    Check whether all subprograms actually exists. Throw otherwise.
    """

    for name, cmd in paths.iteritems():
        if not os.path.isfile( cmd ):
            raise RuntimeError(
                "Subprogram '" + name + "' not found at '" + cmd + "'. Please run setup first."
            )
