# Copyright (C) 2006 Matthias Seeger
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
# USA.

# -------------------------------------------------------------------
# LHOTSE: Toolbox for adaptive statistical models
# -------------------------------------------------------------------
# Module: GLOBAL
# Desc.:  Makefile (generic)
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Customization:
# System-dependent variables are set in the file make.inc.def, or
# alternatively in the file make.inc.$$$ if "make" is called with
# the option "where=$$$". See LHOTSE installation guide for
# details.
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Control variables (pass on command line)
#
# - where: Where does the compilation take place?
#   Use this variable if you maintain the LHOTSE system in different
#   environments in parallel. For a legal value XXX of 'where', a
#   file make.inc.XXX must be present.
#   The default is 'def'. make.inc.def must be present.
# - opt: Level of optimisation:
#		debug:  No optimisation, generate debug info
#		some:   Some optimisation
#		more:   More optimisation. From here on, inlining should
#			be enforced
#		all:	Highest generic level of compiler
#   The default is 'debug'. NOTE: In level 'debug', the compiler does not
#   do inlining, and the matrix/vector methods will run very slow. Do not
#   use this level for final experiments!
# - what: Used only with target 'dep'. Specifies for which main program
#   makedepend should be run (the source file for the main program is
#   main_$(what).cc)
# - mex: Used only for targets resulting in Matlab MEX libraries being
#   built. By def., the value is 'no', and a built will fail when using
#   "mex". This is because LHOTSE enforces static linkage. Values:
#   		yes:	Use this with targets res. in MEX libraries.
#			Will try to use dynamic libraries throughout
#		stat:	Try this if "mex=yes" fails. Here, the static
#			variants of external libraries are used. This
#			requires the variable LOCALLIBS_stat to be set
#			properly in the make.inc.XXX include file
# - prof: Used if code is to be profiled using gprof. By def., the value
#   is 'no'. Profiling code is generated for the value 'yes'
# - fort: By def., the value is 'yes'. If 'no', none of the external
#   Fortran code is compiled. The corr. LHOTSE wrappers will not work,
#   in general an 'NotImplemException' is thrown. The macro HAVE_FORTRAN
#   is defined iff 'fort' is 'yes'
#
# Generic targets:
#
# - depend:    Makedepend for all code
# - dep:       Specific makedepend, using variable 'what' (see above)
# - clean:     Remove all objects, but not the libraries
# - clean_all: Remove all objects and libraries
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# How to include new module:
# - Define xxxDIR,xxxOBJS
# - If module is for main library: add xxxOBJS to OBJECTS
# - Define clean_xxx target, add to clean_int
# - Add entry to DEPENDFILES
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# How to include new main target:
# - define main target xxx (copy-paste from main target list below)
# - define internal main target xxx_int (see below)
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Default values for some of the system-dependent variables
# -------------------------------------------------------------------

CC=		g++
CXX=		g++
EXLIBS_BLAS=	-lf77blas -latlas -lmwlapack -lcblas 
# EXLIBS_BLAS=	-llapack -lf77blas -lcblas -latlas 
FC=		gfortran
# EXLIBS_BLAS=	-llapack -lf77blas -lcblas -latlas -lg2c
# FC=		g77
FFLAGS=		-c -funroll-all-loops -fno-f2c -O3 
MATOPTS=	-cxx
FORCE_STATIC=	-static
MEXCMD=		mex -v -largeArrayDims # -compatibleArrayDims # 

# -------------------------------------------------------------------
# Include make.inc.$(where) file with system-dependent variables
# -------------------------------------------------------------------

mex=	stat
where=	def
prof=	no
opt=	some
fort=	no

include make.inc.$(where)

# -------------------------------------------------------------------
# Optimisation-dependent variables (on 'opt'):
# GCCOPTS:	  Compiler options (machine-indep.). This is the
# 	          concat. of GCCOPTS_gen and GCCOPTS_$(opt), the
#                 latter depends on the optimization level
# LHOTSELIBDIR:   Directory for LHOTSE main library. This is where
#                 the main library is written, and where the source
#                 tree for it is rooted (library modules)
# LHOTSEPROJDIR:  Directory where the LHOTSE projects source tree
#                 is rooted
# LHOTSEEXECDIR:  Directory where LHOTSE executables (of projects)
#                 are written to
# LOCALLIBS:	  External libraries required by LHOTSE
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Compiler options depending on the optimization level in $(opt)
# Inlining is guaranteed to be used from 'more' on. 'more' is a
# good alternative if 'all' requires too much compile time.
# The macro HAVE_INLINE activates inlining with GSL functions.
# -------------------------------------------------------------------

#GCCOPTS_optdebug=	-g -DHAVE_DEBUG
GCCOPTS_optdebug=	-DHAVE_DEBUG -DHAVE_INLINE -O
GCCOPTS_optsome=	-O
GCCOPTS_optmore=	-O2 -funroll-loops -DHAVE_INLINE
GCCOPTS_optall=		-O3 -funroll-loops -DHAVE_INLINE

GCCOPTS_mexno=
GCCOPTS_mexyes=		-DMATLAB_MEX -DUSE_MATLAB_MM -std=c++11 -fPIC 
GCCOPTS_mexstat=	-DMATLAB_MEX -DUSE_MATLAB_MM -std=c++11 -fPIC 

GCCOPTS_profno=
GCCOPTS_profyes=	-pg

LLL_TEMP=	/usr/lib/
LOCALLIBS_mexno=	$(FORCE_STATIC) -lgsl $(EXLIBS_BLAS) -lm
LOCALLIBS_mexyes=	-lgsl $(EXLIBS_BLAS) -lm
LOCALLIBS_mexstat=	$(LLL_TEMP)liblapack.a -lcblas -lf77blas -lgsl -latlas -lgfortran -lm # $(LLL_TEMP)libgsl.a $(LLL_TEMP)libatlas.a # /usr/lib/gcc/x86_64-linux-gnu/4.8/libgfortran.a $(LLL_TEMP)libf77blas.a $(LLL_TEMP)libcblas.a 

INCS_mexno=	-I$(ROOTDIR) $(EXINCS_GSL) $(CXXINCS)
INCS_mexyes=	-I$(ROOTDIR) $(EXINCS_GSL)  $(EXINCS_MATLAB) $(CXXINCS)
INCS_mexstat=	-I$(ROOTDIR) $(EXINCS_GSL)  $(EXINCS_MATLAB) $(CXXINCS)

CXXLDOPTS_profno=
CXXLDOPTS_profyes=	-pg

LHOTSELIBDIR=	$(ROOTDIR)/lhotse
LHOTSEPROJDIR=	$(ROOTDIR)/src
LHOTSEEXECDIR=	$(ROOTDIR)
LHOTSELIB=	$(LHOTSELIBDIR)/liblhotse.a

# -------------------------------------------------------------------
# System-specific options:
# NOT USED IN THE MOMENT!!
# -------------------------------------------------------------------

#LINUXOPTS=      $(GCCOPTS) -march=i686
LINUXOPTS=      $(GCCOPTS)
LINUXLIBS=

SUNOPTS=	$(GCCOPTS)
SUNLIBS=	-lstdc++

ALPHAOPTS=	$(GCCOPTS) -Wno-non-template-friend
ALPHALIBS=	-lstdc++

RANLIB=	ranlib

CPPFLAGS=	$(DEFINES) $(INCS) $(GCCOPTS)

LDFLAGS=	-L$(LHOTSELIBDIR) -L$(LHOTSEPROJDIR) $(EXLDOPTS_BLAS) $(EXLDOPTS_GSL) $(CXXLDOPTS)

LIBS=		$(LOCALLIBS)

# -------------------------------------------------------------------
# Module objects
# -------------------------------------------------------------------

GLOBALDIR=	$(LHOTSELIBDIR)
_GLOBALOBJS=	global \
		StandardException \
		MatlabDebug \
		CommandParser \
		FileUtils \
		MachDep \
		DebugVars \
		IntVal \
		Interval \
		Range \
		ArgBlock \
		TypeCode \
		DefaultLogs \
		FileUtilsMachDep \
		PrintToStdout
#		FixedMemManager
GLOBALOBJS=	$(_GLOBALOBJS:%=$(GLOBALDIR)/%.o)

MATRIXDIR=	$(LHOTSELIBDIR)/matrix
_MATRIXOBJS=	MatDefMembers \
		StVector \
		TempStVector \
		MatStrct \
		StMatrix \
		TempStMatrix \
		IndexVector \
		TempIndexVector \
		LogVector
MATRIXOBJS=	$(_MATRIXOBJS:%=$(MATRIXDIR)/%.o) \

OPTIMIZEDIR=	$(LHOTSELIBDIR)/optimize
_OPTIMIZEOBJS=  OneDimSolver \
		LinearCG \
		GenOpt \
		DoubleLoopOpt \
		ConjGrad \
		QuasiNewton \
		QuasiNewtonDL \
		LimMemBFGSB \
		routines
OPTIMIZEOBJS=	$(_OPTIMIZEOBJS:%=$(OPTIMIZEDIR)/%.o)

QUADDIR=	$(LHOTSELIBDIR)/quad
_QUADOBJS=	GaussQuad \
		GaussProdQuadrature
QUADOBJS=	$(_QUADOBJS:%=$(QUADDIR)/%.o)

RANDODIR=	$(LHOTSELIBDIR)/rando
_RANDOOBJS=	Random
RANDOOBJS=      $(_RANDOOBJS:%=$(RANDODIR)/%.o)

SPECFUNDIR=	$(LHOTSELIBDIR)/specfun
_SPECFUNOBJS=	Specfun \
		InvLogCdfNormal
SPECFUNOBJS=	$(_SPECFUNOBJS:%=$(SPECFUNDIR)/%.o)

MATIFDIR=	$(LHOTSELIBDIR)/matif
_MATIFESSOBJS=	MatlabMatrix \
		MatlabTools
MATIFESSOBJS=	$(_MATIFESSOBJS:%=$(MATIFDIR)/%.o)

ESSENTIALOBJS=	$(GLOBALOBJS) \
		$(MATRIXOBJS) \
		$(RANDOOBJS) \
		$(SPECFUNOBJS) \
		$(OPTIMIZEDIR)/OneDimSolver.o

# -------------------------------------------------------------------
# Project objects
# A project consists of one or more project modules and a main_XXX.cc
# code file. Objects of project modules are not written into the main
# library.
# -------------------------------------------------------------------

EPLINDIR=	$(LHOTSEPROJDIR)/eplin
_EPLINMODOBJS=	ExpectPropLinear \
	     	MaskEPPrior
EPLINMODOBJS=	$(_EPLINMODOBJS:%=$(EPLINDIR)/%.o)

# -------------------------------------------------------------------
# OBJECTS collects all objects for main library LHOTSELIB
# -------------------------------------------------------------------

OBJECTS=	$(GLOBALOBJS) $(ARGBLOCKSOBJS) $(DATAOBJS) \
		$(GPOBJS) $(MATRIXOBJS) $(OPTIMIZEOBJS) \
		$(QUADOBJS) $(RANDOOBJS) $(SPECFUNOBJS)
#		$(EXPEROLDOBJS) $(DATAOLDOBJS) \
#		$(EXPSETUPSOBJS) $(SVMOBJS) $(MYNROBJS) \
#		$(MCMCOBJS) $(DISTOBJS) $(MOFAOBJS) $(SPOTSOBJS) \
#		$(MATLABOBJS) $(INTERFACEOBJS) $(CLASSESOBJS) \
#		$(MLPOBJS)

# -------------------------------------------------------------------
# External Fortran code (used if 'fort' is 'yes'
# -------------------------------------------------------------------

FORTRANOBJS=	$(MATRIXDIR)/dchex.o \
		$(QUADDIR)/gausq2.o

# -------------------------------------------------------------------
# The make process:
#
# 1. The main target <main> defines TARGET as the name of the internal
#     target, which typically is <main>_int, then runs:
# 2. make_opt$(opt), optimisation-specific: This will set the GCCOPTS
#    variable, then TARGET (internal main target).
# 3. make_mex$(mex) (MEX library?): This will set the LOCALLIBS
#    variable, modify the GCCOPTS variable.
#
# NOTE:
# The old Makefile had an automatic recognition of the system architecture:
#
# `uname -s` gives system name, `uname -r` the release.
#
# This is not used anymore.
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Main targets
# 'mainlib' creates the main LHOTSE library, which all projects
# require. The other main targets correspond to projects.
# -------------------------------------------------------------------

epsplinsweep:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

# -------------------------------------------------------------------
# 'opt'-specific  make commands
# 'mex'-specific  make commands
# 'prof'-specific make commands
# -------------------------------------------------------------------

make_optdebug:
	@$(MAKE) make_mex$(mex) GCCOPTS="$(GCCOPTS_optdebug)"

make_optsome:
	@$(MAKE) make_mex$(mex) GCCOPTS="$(GCCOPTS_optsome)"

make_optmore:
	@$(MAKE) make_mex$(mex) GCCOPTS="$(GCCOPTS_optmore)"

make_optall:
	@$(MAKE) make_mex$(mex) GCCOPTS="$(GCCOPTS_optall)"

make_mexno:
	@$(MAKE) make_prof$(prof) GCCOPTS="$(GCCOPTS) $(GCCOPTS_mexno)" LOCALLIBS="$(LOCALLIBS_mexno)" INCS="$(INCS_mexno)"

make_mexyes:
	@$(MAKE) make_prof$(prof) GCCOPTS="$(GCCOPTS) $(GCCOPTS_mexyes)" LOCALLIBS="$(LOCALLIBS_mexyes)" INCS="$(INCS_mexyes)"

make_mexstat:
	@$(MAKE) make_prof$(prof) GCCOPTS="$(GCCOPTS) $(GCCOPTS_mexstat)" LOCALLIBS="$(LOCALLIBS_mexstat)" INCS="$(INCS_mexstat)"

make_profno:
	@$(MAKE) make_fort$(fort) GCCOPTS="$(GCCOPTS) $(GCCOPTS_profno)" CXXLDOPTS="$(CXXLDOPTS_profno)"

make_profyes:
	@$(MAKE) make_fort$(fort) GCCOPTS="$(GCCOPTS) $(GCCOPTS_profyes)" CXXLDOPTS="$(CXXLDOPTS_profyes)"

make_fortno:
	@$(MAKE) $(TARGET)

make_fortyes:
	@$(MAKE) GCCOPTS="$(GCCOPTS) -DHAVE_FORTRAN" OBJECTS="$(OBJECTS) $(FORTRANOBJS)" $(TARGET)

# -------------------------------------------------------------------
# Internal main targets
# -------------------------------------------------------------------

epsplinsweep_int: $(ESSENTIALOBJS) $(QUADOBJS) $(MATIFESSOBJS) $(EPLINMODOBJS)
	$(MEXCMD) $(LHOTSEPROJDIR)/ep_splinsweep.cc $^ $(MATOPTS) CXXFLAGS='$(CPPFLAGS)' $(LDFLAGS) $(LIBS)

# -------------------------------------------------------------------
# Cleaning up stuff
# -------------------------------------------------------------------

clean:
	@$(MAKE) make_opt$(opt) TARGET=$@_int

clean_int:	clean_global clean_matrix clean_optimize \
	clean_quad clean_rando clean_specfun clean_PROJECTMAIN \
	clean_eplin clean_matif

CLEAN_FILES=	-f *.o *.flc *.rpo core errs *~ \#*

clean_LHOTSELIB:
	rm -f $(LHOTSELIB); \

clean_global:
	cd $(GLOBALDIR); \
	rm $(CLEAN_FILES)

clean_matrix:
	@cd $(MATRIXDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_optimize:
	@cd $(OPTIMIZEDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_quad:
	@cd $(QUADDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_rando:
	@cd $(RANDODIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_specfun:
	@cd $(SPECFUNDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_PROJECTMAIN:
	@cd $(LHOTSEPROJDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_eplin:
	@cd $(EPLINDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

clean_matif:
	@cd $(MATIFDIR); \
	rm $(CLEAN_FILES); \
	cd $(ROOTDIR)

# -------------------------------------------------------------------
# Make depend targets
# -------------------------------------------------------------------

DEPENDFILES=	$(GLOBALDIR)/*.cc $(MATRIXDIR)/*.cc \
	$(OPTIMIZEDIR)/*.cc $(QUADDIR)/*.cc $(RANDODIR)/*.cc \
	$(SPECFUNDIR)/*.cc $(MATIFDIR)/*.cc

depend:
	@$(MAKE) make_opt$(opt) TARGET=_depend
#	@$(MAKE) make_`uname -s` REL=`uname -r` TARGET=_depend

dep:
	@$(MAKE) make_opt$(opt) TARGET=_dependspec
#	@$(MAKE) make_`uname -s` REL=`uname -r` TARGET=_dependspec

_depend:
	makedepend -- $(CPPFLAGS) -- $(ROOTDIR)/*.cc $(DEPENDFILES)

_dependspec:
	makedepend -- $(CPPFLAGS) -- $(LHOTSEPROJDIR)/main_$(what).cc \
	$(DEPENDFILES)

# DO NOT DELETE THIS LINE -- make depend depends on it.
