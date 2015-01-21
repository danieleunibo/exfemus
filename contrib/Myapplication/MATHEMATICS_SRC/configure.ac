
AC_INIT(mathematics, 1.0)
AC_CONFIG_AUX_DIR(adm_local)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(mathematics_config.h)

dnl Check Salome Install
CHECK_KERNEL
if test "x$Kernel_ok" = "xno"; then
  AC_MSG_ERROR([You must define a correct KERNEL_ROOT_DIR or use the --with-kernel= configure option !])
fi

dnl Check Salome modules Install


AC_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_CXX
CHECK_F77
CHECK_OMNIORB
CHECK_PACO
CHECK_MPI

MODULE_NAME=mathematics
AC_SUBST(MODULE_NAME)

AC_CHECK_ASTER



echo
echo
echo
echo "------------------------------------------------------------------------"
echo "$PACKAGE $VERSION"
echo "------------------------------------------------------------------------"
echo
echo "Configuration Options Summary:"
echo
echo "  Threads ................ : $threads_ok"
echo "  OmniOrb (CORBA) ........ : $omniORB_ok"
echo "  OmniOrbpy (CORBA) ...... : $omniORBpy_ok"
echo "  Python ................. : $python_ok"
echo "  SALOME KERNEL .......... : $Kernel_ok"
echo "  PaCO++ ................. : $PaCO_ok"
echo "  MPI .................... : $mpi_ok"
echo "  Code Aster ............. : $Aster_ok"

echo
echo "------------------------------------------------------------------------"
echo

if test "x$threads_ok" = "xno"; then
  AC_MSG_ERROR([Thread is required],1)
fi
if test "x$python_ok" = "xno"; then
  AC_MSG_ERROR([Python is required],1)
fi
if test "x$omniORB_ok" = "xno"; then
  AC_MSG_ERROR([OmniOrb is required],1)
fi
if test "x$omniORBpy_ok" = "xno"; then
  AC_MSG_ERROR([OmniOrbpy is required],1)
fi
if test "x$Kernel_ok" = "xno"; then
  AC_MSG_ERROR([SALOME KERNEL is required],1)
fi




AC_CONFIG_FILES([
        Makefile
        idl/Makefile
        resources/Makefile
        src/Makefile
     src/MATH/Makefile
        ])
AC_OUTPUT
