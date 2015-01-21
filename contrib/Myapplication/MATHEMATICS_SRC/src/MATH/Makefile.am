
include $(top_srcdir)/adm_local/make_common_starter.am


lib_LTLIBRARIES= libMATHEngine.la
salomeinclude_HEADERS= MATH_i.hxx


dist_libMATHEngine_la_SOURCES = 	MATH_i.cxx

libMATHEngine_la_CXXFLAGS = -I$(top_builddir)/idl  $(SALOME_INCLUDES) -I${MATHCPP_ROOT_DIR}/include
libMATHEngine_la_LIBADD   = -L${MATHCPP_ROOT_DIR}/lib -lMATHCXX -L$(top_builddir)/idl -lSalomeIDLMATHEMATICS ${SALOME_LIBS} $(FLIBS)



