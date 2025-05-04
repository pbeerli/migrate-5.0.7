AC_INIT([migrate-n], [5.0.6], [beerli@fsu.edu])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign subdir-objects])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL

AC_CHECK_LIB([m], [pow], [MATHLIBS="-lm"], [MATHLIBS=""])
AC_SUBST([MATHLIBS])

AC_ARG_VAR([MPICC], [MPI C compiler])
AC_ARG_VAR([MPICFLAGS], [MPI compiler flags])
AC_ARG_VAR([MPILIBS], [MPI libraries])

if test -n "$MPICC"; then
  AC_MSG_CHECKING([for MPI compile flags])
  MPICFLAGS=`$MPICC -showme:compile`
  MPILIBS=`$MPICC -showme:link`
  echo "+++++++++++++++++++++++++++++++++++++++++++"
  #MPICFLAGS=`$MPICC -show | sed -e 's/^.*gcc//' -e 's/-o.*$$//'`
  AC_MSG_RESULT([$MPICFLAGS])
else
  MPICC=mpicc
  MPICFLAGS=`$MPICC -showme:compile`
  MPILIBS=`$MPICC -showme:link`
  #MPICFLAGS="-I/usr/local/include/mpi"
  #MPILIBS="-lmpi"
fi
AC_SUBST(MPICFLAGS)
AC_SUBST(MPILIBS)

# Always-on preprocessor macros
CPPFLAGS="$CPPFLAGS -DPRETTY -DNEWVERSION -D_REENTRANT"

# Optional: enable debug mode
AC_ARG_ENABLE(debug,
  [AS_HELP_STRING([--enable-debug], [Enable debug build with -g and -Wall])],
  [enable_debug=$enableval],
  [enable_debug=no])
AM_CONDITIONAL([DEBUG_MODE], [test "$enable_debug" = "yes"])
if test "$enable_debug" = "yes"; then
  CFLAGS="-g -Wall -DDEBUG"
  CXXFLAGS="-g -Wall -DDEBUG"
else
  CFLAGS="-O3 -DNDEBUG"
  CXXFLAGS="-O3 -DNDEBUG"
fi
# MPI compiler detection
AC_MSG_CHECKING([for MPI compiler])
AC_CHECK_PROGS([MPICC], [mpicc], [no])
if test "$MPICC" = "no"; then
  AC_MSG_WARN([mpicc not found, MPI build will not work unless you specify MPICC manually])
fi
AC_SUBST(MPICC)

# Check for MPI library (optional but informative)
AC_CHECK_LIB(mpi, MPI_Init, [have_mpi=yes], [have_mpi=no])

# System zlib detection
PKG_CHECK_MODULES([ZLIB], [zlib], [have_zlib=yes], [have_zlib=no])
AM_CONDITIONAL([USE_SYSTEM_ZLIB], [test "$have_zlib" = "yes"])
AC_SUBST([ZLIB_CFLAGS])
AC_SUBST([ZLIB_LIBS])

# System libharu detection
PKG_CHECK_MODULES([LIBHPDF], [libharu], [have_libhpdf=yes], [have_libhpdf=no])
AM_CONDITIONAL([USE_SYSTEM_LIBHPDF], [test "$have_libhpdf" = "yes"])
AC_SUBST([LIBHPDF_CFLAGS])
AC_SUBST([LIBHPDF_LIBS])

# Output
AC_CONFIG_FILES([
  Makefile
  src/Makefile
])
AC_OUTPUT
