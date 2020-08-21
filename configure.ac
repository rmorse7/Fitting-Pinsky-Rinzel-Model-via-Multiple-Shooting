#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)

# Checks for programs.
AC_PROG_CXX
AC_PROG_F77
AC_CHECK_PROG([MEX],[mex],[mex])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lblas':
AC_CHECK_LIB([blas], [main])
# FIXME: Replace `main' with a function in `-lf2c':
AC_CHECK_LIB([f2c], [main])
# FIXME: Replace `main' with a function in `-lginac':
AC_CHECK_LIB([ginac], [main])
echo $LIBS

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stddef.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_C_CONST
#AC_C_INLINE
#AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_ERROR_AT_LINE
#AC_FUNC_MALLOC
#AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([makefile])
AC_CONFIG_FILES([model2diffit])

echo "-- creating symbolic links  --"
cd modules
rm -f *.cc
ln -s ../generic_modules/*.cc .
cd ..
echo "... done."
AC_OUTPUT
chmod uog+x model2diffit

