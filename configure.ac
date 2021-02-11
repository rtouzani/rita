dnl ------------------------------------------------------------------------------
dnl rita - Copyright (C) 2021 Rachid Touzani

dnl This file is part of rita.

dnl rita is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.

dnl rita is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.

dnl ------------------------------------------------------------------------------
dnl Configure script for the rita Software
dnl Process this file with autoconf to produce a configure script.
dnl ------------------------------------------------------------------------------

dnl Disable caching.
define([AC_CACHE_LOAD], )dnl
define([AC_CACHE_SAVE], )dnl

# Define current version and release
sinclude([VERSION])
AC_SUBST(VERSION)
AC_SUBST(RELEASE)

AC_INIT([rita], 1.0.0)
AC_CONFIG_SRCDIR([src/rita.cpp])
AM_INIT_AUTOMAKE

dnl Report what was found
AC_MSG_CHECKING(Version)
AC_MSG_RESULT($VERSION)
AC_MSG_CHECKING(Release)
AC_MSG_RESULT($RELEASE)

AM_CONFIG_HEADER(src/config.h)
AM_MAINTAINER_MODE

AC_CYGWIN()

dnl Prevents setting flags.
CPPFLAGS="-O3 -Wall "
CXXFLAGS=" "
CFLAGS="-O3 -Wall "

if test "$OFELI_CFLAGS" == "" ; then
   OFELI_CFLAGS="-I$includedir/ofeli "
fi
if test "$OFELI_LIBS" == "" ; then
   OFELI_LIBS="-lofeli "
fi

if test "$GMSH_CFLAGS" == "" ; then
   GMSH_CFLAGS="-I$includedir "
fi
if test "$GMSH_LIBS" == "" ; then
   GMSH_LIBS="-lgmsh "
fi

AC_SUBST(RITA_CFLAGS)
AC_SUBST(OFELI_CFLAGS)
AC_SUBST(OFELI_LIBS)
AC_SUBST(OFELI_LDFLAGS)
AC_SUBST(GMSH_CFLAGS)
AC_SUBST(GMSH_LIBS)
AC_SUBST(GMSH_LDFLAGS)
AC_SUBST(SYMBOLIC_CFLAGS)


dnl Checks for programs.
AC_LANG_CPLUSPLUS
AC_PROG_MAKE_SET
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL

dnl Get the operating system and machine names
UNAME=`uname`
HOSTNAME=`hostname`

dnl Check compiler features
AC_TYPE_SIZE_T
AC_C_CONST
AC_C_INLINE

dnl Add copious amounts of debugging with gcc, egcs, etc
if test "${GXX}" = "yes" ; then
CPPFLAGS="${CPPFLAGS}"
fi
if test "${GCC}" = "yes" ; then
CFLAGS="${CFLAGS}"
fi

AC_ARG_VAR([RITA_CFLAGS], [rita Compilation flags])
AC_ARG_VAR([OFELI_CFLAGS], [OFELI Compilation flags.
                            Default = "-I${includedir}/ofeli"])
AC_ARG_VAR([OFELI_LIBS], [OFELI Library option.
                         Default = "-lofeli"])
AC_ARG_VAR([GMSH_CFLAGS], [Gmsh Compilation flags.
                           Default = "-I${includedir}"])
AC_ARG_VAR([GMSH_LIBS], [Gmsh Library option.
                         Default = "-lgmsh"])
AC_ARG_VAR([SYMBOLIC_CFLAGS], [Symbolicc++ Compilation flags.
                               Default = "-Isymbolicc++"])

CPPFLAGS="${CPPFLAGS} ${OFELI_CFLAGS} ${SYMBOLIC_CFLAGS}"
LIBS="$LIBS $OFELI_LIBS $OFELI_LDFLAGS"

dnl Building with GMSH
AC_MSG_CHECKING(for Gmsh)
AC_ARG_WITH(gmsh,[  --with-gmsh             use Gmsh])
AC_MSG_RESULT([$with_gmsh])
if test "x$with_gmsh" = "xyes" ; then
   CPPFLAGS="${CPPFLAGS} -DUSE_GMSH ${GMSH_CFLAGS}"
   LIBS="$LIBS -DUSE_GMSH $GMSH_LIBS $GMSH_LDFLAGS"
#   AC_CHECK_LIB(gmsh, send, [GMSH_LIBS], [])
fi
dnl Debugging turned on
AC_MSG_CHECKING(for debugging)
AC_ARG_ENABLE(debug,[  --enable-debug          compile for debugging])
AC_MSG_RESULT([$enable_debug])
if test "x$enable_debug" = "xyes" ; then
CPPFLAGS="${CPPFLAGS} -g -DDEBUG"
CFLAGS="${CFLAGS} -g -DDEBUG"
fi

dnl Debugging and Tracing turned on
AC_MSG_CHECKING(for debugging and tracing)
AC_ARG_ENABLE(debug,[  --enable-trace-debug    compile for debugging and tracing])
AC_MSG_RESULT([$enable_trace_debug])
CPPFLAGS="${CPPFLAGS}"
CFLAGS="${CFLAGS}"
if test "x$enable_trace_debug" = "xyes" ; then
CPPFLAGS="${CPPFLAGS} -g -DDEBUG "
CFLAGS="${CFLAGS} -g -DDEBUG "
fi

dnl Building for release
AC_MSG_CHECKING(for release build)
AC_ARG_ENABLE(release,[  --enable-release        compile for release])
AC_MSG_RESULT([$enable_release])
if test "x$enable_release" = "xyes" ; then
CPPFLAGS="${CPPFLAGS} -O3 -DNDEBUG"
CFLAGS="${CFLAGS} -O3 -DNDEBUG"
if test "${GXX}" = "yes" ; then
CPPFLAGS="${CPPFLAGS} -ffast-math -fstrict-aliasing"
CFLAGS="${CFLAGS} -ffast-math -fstrict-aliasing"
fi
fi
 
dnl Building for fast
AC_MSG_CHECKING(for fast build)
AC_ARG_ENABLE(fast,[  --enable-fast           compile for fast configuration])
AC_MSG_RESULT([$enable_fast])
if test "x$enable_fast" = "xyes" ; then
CPPFLAGS="${CPPFLAGS} -Ofast -DNDEBUG"
CFLAGS="${CFLAGS} -Ofast -DNDEBUG"
if test "${GXX}" = "yes" ; then
CPPFLAGS="${CPPFLAGS} -ffast-math -fstrict-aliasing"
CFLAGS="${CFLAGS} -ffast-math -fstrict-aliasing"
fi
fi

# Restore CFLAGS
CFLAGS="$save_cflags "

# Variables that will be substituted in the Makefile.in's

AC_SUBST(PACKAGE)
inc_dir="$includedir"

# Now restore our (possibly augmented) CFLAGS, LIBS and CPPFLAGS.
CFLAGS="$save_cflags"
RITA_CFLAGS="$RITA_CFLAGS"
RITA_LIBS="$OFELI_LIBS $OFELI_LDFLAGS"
if test "x$with_gmsh" = "xyes" ; then
RITA_LIBS="$GMSH_LIBS $GMSH_LDFLAGS"
fi

AC_OUTPUT(Makefile
          src/Makefile
          tutorial/Makefile
          tutorial/ae/Makefile
          tutorial/ode/Makefile
          tutorial/pde/Makefile
          tutorial/eigen/Makefile
          tutorial/optim/Makefile
          tutorial/integration/Makefile
          tutorial/approximation/Makefile
          doc/Makefile
          doc/css/Makefile
          doc/fonts/Makefile
          doc/fonts/bootstrap/Makefile
          doc/fonts/icomoon/Makefile
          doc/images/Makefile
          doc/js/Makefile
          doc/sass/Makefile
          doc/sass/bootstrap/Makefile
          doc/sass/bootstrap/mixins/Makefile
          test_rita.sh
)

AC_SUBST(PACKAGE)

echo "********************************************************************"
echo "rita is configured for"
echo "  - OS                  : ${UNAME} on ${HOSTNAME}"
echo "  - C++ compiler        : ${CXX}"
echo "  - Exe linker          : ${CXX}"
echo "  - C++ flags           : ${CPPFLAGS}"
echo "********************************************************************"
echo ""
echo "Configure finished!"
echo "            Type: 'make' to compile rita."
echo "            Then: 'make install' to install it."
echo "            Type: 'make tests' to run tutorial tests"
echo ""
