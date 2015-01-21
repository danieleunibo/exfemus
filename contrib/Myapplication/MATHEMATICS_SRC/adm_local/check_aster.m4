
#
# Check availability of Aster binary distribution
#

AC_DEFUN([AC_CHECK_ASTER],[

AC_CHECKING(for Aster)

Aster_ok=no

AC_ARG_WITH(aster,
      [AC_HELP_STRING([--with-aster=DIR],[root directory path of Aster installation])],
      [ASTER_DIR="$withval"],[ASTER_DIR=""])

if test -f ${ASTER_DIR}/bin/aster ; then
   Aster_ok=yes
   AC_MSG_RESULT(Using Aster distribution in ${ASTER_DIR})

   ASTER_INCLUDES=-I$ASTER_DIR/include/aster

   AC_SUBST(ASTER_DIR)
   AC_SUBST(ASTER_INCLUDES)

else
   AC_MSG_WARN("Cannot find Aster distribution")
fi

AC_MSG_RESULT(for Aster: $Aster_ok)

])dnl
