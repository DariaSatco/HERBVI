# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Jimmy],[4.31],[jimmy@cedar.ac.uk],[jimmy])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE
#AC_LANG(Fortran 77)
AC_LANG(C++)

AC_CONFIG_SRCDIR([src/divon4/addbnd.F])
AC_CONFIG_HEADERS([config/config.h])

## Checks for programs.
AC_PROG_F77
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_LN_S

# Check for required Herwig
AC_CEDAR_CHECK_GENHEAD([Herwig], [6510])

AC_CONFIG_FILES([Makefile src/Makefile include/Makefile])
AC_OUTPUT

