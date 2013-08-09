
AC_PREREQ(2.59)
AC_INIT([CoalSweep], [0.01], [mtholder@ku.edu])

AC_CONFIG_SRCDIR([src/chimne_ssweep_main.cpp])
AC_CONFIG_HEADER([config.h])

#	Directory that contains install-sh and other auxiliary files
#AC_CONFIG_AUX_DIR([config])

################################################################################
# Checks for build-platform and target info
#	this defines the "target" variable that is used later in this file
################################################################################
AC_CANONICAL_TARGET

################################################################################
#	According to (http://www.mail-archive.com/autoconf@gnu.org/msg14232.html)
#		this macro should be after AC_INIT but before AM_INIT_AUTOMAKE
################################################################################
#AC_CONFIG_MACRO_DIR(config/m4)

AM_INIT_AUTOMAKE([1.9 foreign dist-tarZ tar-ustar filename-length-max=299])

# Checks for programs.
INP_CFLAGS="$CFLAGS"
AC_PROG_CC
AC_PROG_CPP
INP_CXXFLAGS="$CXXFLAGS"
AC_PROG_CXX

# Initialize CXXFLAGS to prevent it from defaulting to "-g -O2"
CXXFLAGS="$INP_CXXFLAGS"
CFLAGS="$CFLAGS $INP_CFLAGS"

# required because we are building a library 
#AC_PROG_RANLIB


AC_PROG_INSTALL
#AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h malloc.h stddef.h stdlib.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr strdup strtol])



#--------------------------------------------------------------------------------------------------#
#  Set CXXFLAGS                                                                                    #
#--------------------------------------------------------------------------------------------------#


# Check for debugging mode.
AC_ARG_ENABLE(debugging,
	AC_HELP_STRING([--enable-debugging],[build for debugging]), , [enable_debugging=no])
if test "$enable_debugging" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: debugging is enabled; optimization is suppressed!
	   ])
fi

# Check whether asserts should be allowed.
AC_ARG_ENABLE(asserts,
	AC_HELP_STRING([--enable-asserts],[build with asserts on (NDEBUG not defined)]), , [enable_asserts=no])
if test "$enable_asserts" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: compiling with assertions on (NDEBUG not defined)
	   ])
fi

# Compile openMP multithreaded version
AC_ARG_ENABLE(openmp,
	AC_HELP_STRING([--enable-openmp],[build OpenMP multithreaded version]), , [enable_openmp=no])
if test "$enable_openmp" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: compiling OpenMP multithreaded version
	   ])
fi

# Mpi run forking version
AC_ARG_ENABLE(mpi,
	AC_HELP_STRING([--enable-mpi],[build MPI run distributing version]), , [enable_mpi=no])
if test "$enable_mpi" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: compiling MPI run distributing version
	   ])
fi

# Initialize optimization flag in case it doesn't get set below.
CXXFLAGS_OPTIM_SPEED="-O"

# "-g" may not work with some compilers, but end users shouldn't be
if test "$enable_debugging" = yes; then
	CXXFLAGS_OPTIM_SPEED="-O0"
	CXXFLAGS="$CXXFLAGS -Wall -g -Wreturn-type -Wunused -Wredundant-decls -Wcast-align -Wcomment -Wextra"
fi

#DJZ allow turning off NDEBUG only
if test  "$enable_asserts" = no; then
	CXXFLAGS="$CXXFLAGS -DNDEBUG"
fi

#DJZ OpenMP compilation
if test  "$enable_openmp" = yes; then
	if test "$CC" = "icc"; then
		CXXFLAGS="$CXXFLAGS -openmp"
	else
		AC_MSG_ERROR([the --enable-openmp option can only be used with the intel compiler (CC=icc)])
	fi
fi

if test "$CC" = "icc"; then
#	Intel C compiler for Linux
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O2 -ip -funroll-loops -fno-alias"
	fi
elif test "$CC" = "ccc"; then
#	Compaq C compiler for Linux
	if test "x$arch" = "x"; then
		arch="host"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-fast -inline speed -arch $arch"
	fi
elif test "$CC" = "xlc"; then
#	IBM XL C compiler
	CFLAGS="$CXXFLAGS -qsourcetype=c++ -qenablevmx -qchars=signed"
	if test "x$arch" = "x"; then
		arch="auto"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O3 -qarch=$arch -qtune=$arch -qalias=ansi -qunroll=yes"
	fi
elif test "x$GCC" = "xyes" ; then
	CXXFLAGS="$CXXFLAGS -Wno-uninitialized"
	if test "$enable_debugging" = yes; then
		CXXFLAGS_OPTIM_SPEED="-O0 -Wimplicit"
	else
		CXXFLAGS_OPTIM_SPEED="-O3 -ffast-math -funroll-loops -fstrict-aliasing"
	fi
	case "$build_os" in
	 	darwin*) CXXFLAGS="$CXXFLAGS -Wno-long-double" ;;
	 	*) CXXFLAGS="$CXXFLAGS -fsigned-char";;
	esac	
fi

AC_LANG(C++)
ACX_MPI

#MPI run forker
if test  "$enable_mpi" = yes; then
	CXXFLAGS="$CXXFLAGS -DSUBROUTINE_COALSWEEP"
	CC="$MPICC"
	CXX="$MPICXX"
	LIBS="$MPILIBS $LIBS"
fi
CXXFLAGS="$CXXFLAGS $CXXFLAGS_OPTIM_SPEED"


################################################################################
#	By default, use the bundled ncl code from source, 
#		but let the --with-ncl argument override this.
################################################################################
NCL_INC_DIR="/usr/local/include"
NCL_LIB_DIR="/usr/local/lib"
AC_ARG_WITH(
	[ncl],
	AC_HELP_STRING(
		[--with-ncl=DIR],
		[Specify the root directory for ncl library (parent of the include/ncl and lib directories). Omit this flag to build and use the bundled version of ncl.]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		NCL_INC_DIR="$withval/include"
		NCL_LIB_DIR="$withval/lib"
	fi
	])

if ! test -d "$NCL_INC_DIR/ncl" ; then
	AC_MSG_ERROR([NCL is a prerequisite for building coalsweep. Expecting a directory called $NCL_INC_DIR/ncl, verify that the argument to --with-ncl is correct. You may need to manually make the directory $NCL_INC_DIR/ncl copy all of the ncl headers into that directory if you are using NCL 2.0])
else
	if ! test -f "$NCL_INC_DIR/ncl/nxsdefs.h" ; then
		AC_MSG_ERROR([NCL is a prerequisite for building coalsweep. Expecting a directory called $NCL_INC_DIR/ncl, verify that the argument to --with-ncl is correct. You may need to manually make the directory $NCL_INC_DIR/ncl copy all of the ncl headers into that directory if you are using NCL 2.0])
	fi
fi

if ! test -d "$NCL_LIB_DIR" ; then
	AC_MSG_ERROR([NCL is a prerequisite for building coalsweep. Expecting a directory called $NCL_LIB_DIR, verify that the argument to --with-ncl is correct.])
fi

INCLUDE="$INCLUDE -I$NCL_INC_DIR -I$NCL_INC_DIR/ncl"
if test -d "$NCL_LIB_DIR/ncl" ; then
	LDFLAGS="$LDFLAGS -L$NCL_LIB_DIR/ncl"
else
	LDFLAGS="$LDFLAGS -L$NCL_LIB_DIR"
fi


CHIMNE_SSWEEP_LIBS="-lncl"

C_PHY_PROB_DIR=""
AC_ARG_WITH(
	[cphyprob],
	AC_HELP_STRING(
		[--with-cphyprob=DIR],
		[Specify path to the  parent of the cPhyProb/ccore directory which contains dsct_model.h]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		C_PHY_PROB_DIR="$withval"
	fi
	])

if ! test -f "$C_PHY_PROB_DIR/cPhyProb/ccore/dsct_model.h" ; then
	AC_MSG_ERROR([cPhyProb is a prerequisite for building chimnesweep. Expecting a file at  $C_PHY_PROB_DIR/cPhyProb/ccore/dsct_model.h, verify that the argument to --with-cphyprob is correct.  ])
fi


BEAGLE_DIR=""
AC_ARG_WITH(
	[beagle],
	AC_HELP_STRING(
		[--with-beagle=DIR],
		[Specify path to the  parent of the beagle prefix directory which contains the include and lib directories]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		BEAGLE_DIR="$withval"
	fi
	])

if ! test -z "$BEAGLE_DIR"
then
	if ! test -f "$BEAGLE_DIR/include/libhmsbeagle-1.0/libhmsbeagle/beagle.h" ; then
		AC_MSG_ERROR([Expecting a file at $BEAGLE_DIR/include/libhmsbeagle-1.0/libhmsbeagle/beagle.h, verify that the argument to --with-beagle is correct.  ])
	fi
	CFLAGS="$CFLAGS -DUSE_BEAGLE_LIB"
	CPPFLAGS="$CPPFLAGS -I$BEAGLE_DIR/include/libhmsbeagle-1.0"
	CXXFLAGS="$CXXFLAGS -DUSE_BEAGLE_LIB"
	LDFLAGS="$LDFLAGS -L$BEAGLE_DIR/lib"
	CHIMNE_SSWEEP_LIBS="$CHIMNE_SSWEEP_LIBS -lhmsbeagle"
fi


CPPFLAGS="$CPPFLAGS -I$NCL_INC_DIR -I$NCL_INC_DIR/ncl -I$C_PHY_PROB_DIR -DEDGE_LEN_FROM_DEPTH"
if test -d "$NCL_LIB_DIR/ncl" ; then
	LDFLAGS="$LDFLAGS -L$NCL_LIB_DIR/ncl"
else
	LDFLAGS="$LDFLAGS -L$NCL_LIB_DIR"
fi


AC_SUBST([CHIMNE_SSWEEP_LIBS])

AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])


#currently we are bundling NCL in the svn, so we know the path

# Variables to substitute in Makefile:
AC_SUBST(CPPFLAGS)

AC_CONFIG_FILES([
	Makefile 			\
	src/Makefile		\
	])

AC_OUTPUT
