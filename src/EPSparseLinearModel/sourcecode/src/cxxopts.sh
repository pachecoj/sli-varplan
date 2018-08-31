#
# cxxopts.sh	Shell script for configuring MEX-file creation script,
#               mex.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  The version tested against
#               is specified by platform.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# Copyright 1984-2000 The MathWorks, Inc.
# $Revision: 1.51 $  $Date: 2000/06/13 14:41:45 $
#----------------------------------------------------------------------------
#
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The cmex script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        alpha)
#----------------------------------------------------------------------------
            #
            # DIGITAL C++ V6.0-010
            #
            CC='cxx'
            CFLAGS=''
            CLIBS=''
            COPTIMFLAGS='-O2 -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='-shared'
            FLIBS='-lUfor -lfor -lFutil'
            FOPTIMFLAGS='-O2'
            FDEBUGFLAGS='-g'
#
            LD="cxx"
            LDFLAGS="-shared -Wl,-expect_unresolved,'*',-hidden,-exported_symbol,'__*',-exported_symbol,$ENTRYPOINT,-exported_symbol,mexVersion"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        hpux)
#----------------------------------------------------------------------------
            #
            # HP aC++ B3910B A.01.06
            #
            CC='aCC'
            CFLAGS='+z -D_POSIX_C_SOURCE=199506L -D_XOPEN_SOURCE'
            CLIBS=""
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f90'
            FFLAGS='+z'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='aCC'
            LDFLAGS="-b -Wl,+e,$ENTRYPOINT,+e,mexVersion,+e,_shlInit"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        hp700)
#----------------------------------------------------------------------------
            #
            # HP aC++ B3910B A.01.06
            #
            CC='aCC'
#
# +DAportable - remove from CFLAGS if you wish to optimize for target machine
#
            CFLAGS='+z -D_HPUX_SOURCE +DAportable'
            CLIBS=''
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f90'
#
# +DAportable - remove from FFLAGS if you wish to optimize for target machine
#
            FFLAGS='+z +DAportable'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='aCC'
            LDFLAGS="-b -Wl,+e,$ENTRYPOINT,+e,mexVersion,+e,_shlInit"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        ibm_rs)
#----------------------------------------------------------------------------
            #
            # AIX Compilers V3.1
            #
            CC='xlC'
            CFLAGS='-DNO_BUILT_IN_SUPPORT_FOR_BOOL'
            CLIBS="-L$MATLAB/bin/$Arch -lmex -lmx -lmatlbmx -lm -lC"
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS=''
            FLIBS="$MATLAB/extern/lib/ibm_rs/fmex1.o -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='xlC'
            LDFLAGS="-bI:$MATLAB/extern/lib/ibm_rs/exp.ibm_rs -bE:$MATLAB/extern/lib/ibm_rs/$MAPFILE -bM:SRE -e $ENTRYPOINT"
            LDOPTIMFLAGS='-s'
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        glnx86)
#----------------------------------------------------------------------------
            CC='g++'
#
# Exception handling does not work with g++ version, V2.7.2.1
#
#           CFLAGS='-fPIC -fhandle-exceptions'
#
            CFLAGS='-fPIC'
            CLIBS=''
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='g77'
            FFLAGS='-fPIC'
            FLIBS='-lf2c'
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD=$CC
            LDFLAGS='-shared'
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        sgi)
#----------------------------------------------------------------------------
            #
            # Base Compiler Development Environment, 7.2.1
            #
            CC='CC'
            #
            # Add -exceptions to CFLAGS to support exception handling
            #
            CFLAGS='-n32 -mips3'
            CLIBS="-lC"
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='-n32 -mips3'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='CC'
            LDFLAGS="-shared -exported_symbol $ENTRYPOINT -exported_symbol mexVersion"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        sol2)
#----------------------------------------------------------------------------
            #
            # WorkShop Compilers 5.0
            #
            CC='CC'
            CCV=`CC -V 2>&1`
            version=`expr "$CCV" : '.*\([0-9][0-9]*\)\.'`
            if [ "$version" = "4" ]; then
                    echo "SC5.0 or later C++ compiler is required"
            fi
            CFLAGS='-KPIC -dalign'
            CLIBS='-lCstd -lCrun'
	    COPTIMFLAGS='-O -DNDEBUG'
	    CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='-dalign'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD=$CC
            LDFLAGS="-G -M $MATLAB/extern/lib/sol2/$MAPFILE"
	    LDOPTIMFLAGS=''
	    LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
