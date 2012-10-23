## VBmix
AC_INIT([VBmix], [0.2.10], [pbruneau@gmail.com])
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi


CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`




AC_CHECK_LIB(m, cos, , AC_MSG_ERROR("math library not found"))
AC_CHECK_LIB(gslcblas, cblas_dgemm, , 
	AC_MSG_ERROR("gslcblas library not found"))
AC_CHECK_LIB(gsl, gsl_blas_dgemm, , 
	AC_MSG_ERROR("gsl library not found"))
AC_CHECK_LIB(fftw3, fftw_execute, , 
	AC_MSG_ERROR("fftw3 library not found"))
AC_CHECK_LIB(QtCore, main, [HAVE_QT=1], [HAVE_QT=0])
if test ${HAVE_QT} = 0; then
	AC_CHECK_HEADER(/Library/Frameworks/QtCore.framework/Headers/qconfig.h, , AC_MSG_ERROR("QtCore library not found"))
fi

: ${TP_ARCH=`uname -s`}
case "${TP_ARCH}" in
	*Darwin*)
	CPPFLAGS="$CPPFLAGS -I/Library/Frameworks/QtCore.framework/Headers"
	LIBS="$LIBS -framework QtCore -lfftw3 -lgsl -lgslcblas -lm"
	;;
	*)
	QTFLAGS=`pkg-config QtCore --cflags`
	CPPFLAGS="$CPPFLAGS $QTFLAGS"
	LIBS="$LIBS -lQtCore -lfftw3 -lgsl -lgslcblas"
	;;
esac

#  	*Linux*)
#	CPPFLAGS="$CPPFLAGS -I/usr/include/qt4 -I/usr/include/qt4/QtCore"
#	LIBS="$LIBS -lQtCore -lfftw3 -lgsl -lgslcblas"
#	;;
#  *)
#	AC_MSG_ERROR([Unsupported operating system: ${TP_ARCH}])
#	;;



AC_SUBST(LIBS)
AC_SUBST(CPPFLAGS)
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

