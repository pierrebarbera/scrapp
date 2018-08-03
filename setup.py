#!/usr/bin/python

import os
import scripts.util as util

def fatal(msg=""):
    if msg:
        print msg
    print "Terminating..."
    exit(1)

# create deps/bin folder

# detect existence of installed programs
# git
# flex
# bison
#


# TODO ensure setups are always rerun completely, incase the user updates the submodule

print "Detected machine characteristics:"
mach = util.get_platform()
print mach

# get all subprograms
progs = util.subprogram_commands()

# iterate over them, check if they are present, if not, fetch them
for name, cmd in progs.iteritems():
    print "Checking for " + name + "..."
    if not os.path.isfile( cmd ):
        print name + " was not found at path " + cmd + "!"
        if ( util.try_resolve(name) != 0 ):
            raise RuntimeError( "Failed to resolve " + name )
