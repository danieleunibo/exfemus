dnl modify CFLAGS, CXXFLAGS and LIBS for compiling pthread-based programs.
dnl@author  (C) Ruslan Shevchenko <Ruslan@Shevchenko.Kiev.UA>, 1998, 2000
dnl@id  $Id: enable_pthreads.m4,v 1.7 2012-08-08 13:47:08 vsr Exp $
dnl Modified to use acx_pthread.m4 from GNU Autoconf Macro Archive
dnl
AC_DEFUN([ENABLE_PTHREADS],[
AC_REQUIRE([ACX_PTHREAD])

if test x"$enable_pthreads_done" != xyes; then
  if test x"$acx_pthread_ok" = xyes; then
    CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"
    LIBS="$LIBS $PTHREAD_LIBS"
    threads_ok=yes
  else
    threads_ok=no
  fi
  enable_pthreads_done=yes
fi
])dnl
dnl
