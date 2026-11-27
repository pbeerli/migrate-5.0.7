AC_INIT([migrate-n], [6.1.112], [beerli@fsu.edu])
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

# libharu: the bundled lib/haru (synced from libharu-pdfua, which carries
# fixes that tagged PDF reports need) unless --with-system-libhpdf. A system
# libharu is tested for the /ToUnicode CMap fix (tools/check_libhpdf_tounicode.c);
# --with-system-libhpdf=force links it without that test.
AC_ARG_WITH([system-libhpdf],
  [AS_HELP_STRING([--with-system-libhpdf@<:@=force@:>@],
    [link a system libharu (pkg-config libhpdf) instead of the bundled one; it is
     checked for the /ToUnicode fix of libharu-pdfua, force skips the check @<:@no@:>@])],
  [], [with_system_libhpdf=no])
have_libhpdf=no
if test "x$with_system_libhpdf" != xno; then
  PKG_CHECK_MODULES([LIBHPDF], [libhpdf], [have_libhpdf=yes],
    [AC_MSG_ERROR([--with-system-libhpdf: pkg-config finds no libhpdf; drop the option to use the bundled libharu])])
  if test "x$with_system_libhpdf" = xforce; then
    AC_MSG_WARN([system libharu linked without the /ToUnicode check (--with-system-libhpdf=force)])
  elif test "x$cross_compiling" = xyes; then
    AC_MSG_WARN([cross-compiling: cannot run the /ToUnicode check of the system libharu])
  else
    AC_MSG_CHECKING([whether the system libharu writes a valid /ToUnicode CMap])
    hpdf_check_out=`$CC $CFLAGS $LIBHPDF_CFLAGS -I$srcdir/src -o conftest_hpdf$EXEEXT $srcdir/tools/check_libhpdf_tounicode.c $LDFLAGS $LIBHPDF_LIBS -lz -lm 2>&1 && ./conftest_hpdf$EXEEXT 2>&1`
    hpdf_check_rc=$?
    rm -f conftest_hpdf$EXEEXT
    case $hpdf_check_rc in
      0) AC_MSG_RESULT([yes]) ;;
      1) AC_MSG_RESULT([no])
         AC_MSG_ERROR([the system libharu lacks the /ToUnicode CMap fix that libharu-pdfua
  carries ($hpdf_check_out). Tagged PDF reports would not map symbol glyphs
  to Unicode and would fail PDF/UA-1. Drop --with-system-libhpdf to use the
  bundled libharu, or use --with-system-libhpdf=force to link it anyway.]) ;;
      *) AC_MSG_RESULT([unknown])
         AC_MSG_ERROR([could not build or run tools/check_libhpdf_tounicode.c against the
  system libharu: $hpdf_check_out
  (--with-system-libhpdf=force links it without the check)]) ;;
    esac
  fi
fi
AM_CONDITIONAL([USE_SYSTEM_LIBHPDF], [test "$have_libhpdf" = "yes"])
AC_SUBST([LIBHPDF_CFLAGS])
AC_SUBST([LIBHPDF_LIBS])
if test "$have_libhpdf" = yes; then
  AC_MSG_NOTICE([libharu: system library ($LIBHPDF_LIBS)])
else
  hpdf_synced=`sed -n 's/.*libharu-pdfua commit \([[0-9a-f]]*\).*/libharu-pdfua \1/p' $srcdir/lib/haru/SYNCED_FROM 2>/dev/null`
  AC_MSG_NOTICE([libharu: bundled lib/haru ($hpdf_synced)])
fi

# Output
AC_CONFIG_FILES([
  Makefile
  src/Makefile
  lib/haru/include/hpdf_config.h
])
AC_OUTPUT

