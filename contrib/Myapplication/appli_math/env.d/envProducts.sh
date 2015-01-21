##########################################################################
# Functions exporta and exportp are used to append/prepend correspondingly
# one directory or a set of directories separated by semicolon symbol (':')
# to the environment variables like PATH, LD_LIBRARY_PATH, PYTHONPATH,
# LD_RUN_PATH etc.
# The main purpose is to replace default setenv command behavior:
# exporta, exportp also remove duplicated entries, shortening in that way
# the environment variables.
# If some directory being added is already included into the variable
# nothing is done for it.
# Note, that these functions work some slower that setenv command itself.
#
#### cleandup ###
# appends/prepends set of directories (second parameter)
# to the another set of directories (first parameter) and
# removes duplicated entries;
# the third parameter defines the mode: 0 - append, 1 - prepend
cleandup() {
out_var=`echo $1 $2 | awk -v dir=$3 '{                   \
     na = split($2,a,":");                               \
     k1=0;                                               \
     bbb[k1++]="";                                       \
     ccc[""];                                            \
     if($1 != "<empty>") {                               \
       nb = split($1,b,":");                             \
       for(i=1;i<=nb;i++) {                              \
         if(!(b[i] in ccc) ) {                           \
	   ccc[b[i]];                                    \
           bbb[k1++]=b[i];                               \
	 };                                              \
       };                                                \
     };                                                  \
     k2=0;                                               \
     aaa[k2++]="";                                       \
     for(i=1;i<=na;i++) {                                \
       if(!(a[i] in ccc)) {                              \
         ccc[a[i]];                                      \
         aaa[k2++]=a[i];                                 \
       };                                                \
     };                                                  \
     ORS=":";                                            \
     if(dir) {                                           \
       for(i=1;i<k2;i++) {                               \
         print aaa[i];                                   \
       }                                                 \
       for(i=1;i<k1;i++) {                               \
         print bbb[i];                                   \
       }                                                 \
     }                                                   \
     else {                                              \
       for(i=1;i<k1;i++) {                               \
         print bbb[i];                                   \
       }                                                 \
       for(i=1;i<k2;i++) {                               \
         print aaa[i];                                   \
       }                                                 \
     }                                                   \
   }' | sed -e 's/\(.*\):/\1/g'`
echo ${out_var}
}
### exporta ###
# appends directory or set of directories, separated by ':' (second parameter)
# to the variable (first parameter)
exporta () {
xenv=${!1}
if [ -z "${xenv}" ]; then xenv="<empty>"; fi
out_var=`cleandup ${xenv} $2 0`
export $1=${out_var}
}
### exportp ###
# prepends directory or set of directories, separated by ':' (second parameter)
# to the variable (first parameter)
exportp () {
xenv=${!1}
if [ -z "${xenv}" ]; then xenv="<empty>"; fi
out_var=`cleandup ${xenv} $2 1`
export $1=${out_var}
}
###########################################################################

#------ Setting products installation directory ------
export INST_ROOT=/homesd/msandro/software/salome

#------ Environment switch: 0 for SALOME building, 1 for SALOME launching ------
export ENV_FOR_LAUNCH=1

#------ Setting numeric locale ------
export LC_NUMERIC=C

#------ tcltk ------
export TCLHOME=${INST_ROOT}/tcltk-8.6.0
export TCLLIBPATH="${TCLHOME}/lib ${TCLHOME}/lib/tcl8.6 ${TCLHOME}/lib/tk8.6"
exportp PATH ${TCLHOME}/bin
exportp LD_LIBRARY_PATH ${TCLHOME}/lib
##
## #------ tcltk_src ------
## # nothing to do
## ##
#------ Python ------
export PYTHON_ROOT_DIR=${INST_ROOT}/Python-2.7.3
export PYTHONHOME=${PYTHON_ROOT_DIR}
export PYTHON_VERSION=2.7
export PYTHON_INCLUDE=${PYTHON_ROOT_DIR}/include/python${PYTHON_VERSION}
exportp PATH ${PYTHON_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${PYTHON_ROOT_DIR}/lib
exportp PYTHONPATH ${PYTHON_ROOT_DIR}/lib/python${PYTHON_VERSION}
##
## #------ Python_src ------
## # nothing to do
## ##
#------ Qt ------
export QT4_ROOT_DIR=${INST_ROOT}/qt-4.8.4
export QT_PLUGIN_PATH=${QT4_ROOT_DIR}/plugins
exportp PATH ${QT4_ROOT_DIR}/bin 
exportp LD_LIBRARY_PATH ${QT4_ROOT_DIR}/lib
##
## #------ Qt_src ------
## # nothing to do
## ##
#------ Sip ------
export SIP_ROOT_DIR=${INST_ROOT}/sip-4.14.2
exportp PATH ${SIP_ROOT_DIR}/bin
exportp PYTHONPATH ${SIP_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages
exportp LD_LIBRARY_PATH ${SIP_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages
##
## #------ Sip_src ------
## # nothing to do
## ##
#------ PyQt ------
export PYQT4_ROOT_DIR=${INST_ROOT}/PyQt-4.9.6
export PYQT_SIPS=${PYQT4_ROOT_DIR}/sip
exportp PATH ${PYQT4_ROOT_DIR}/bin
exportp PYTHONPATH ${PYQT4_ROOT_DIR}:${PYQT4_ROOT_DIR}/PyQt4
exportp LD_LIBRARY_PATH ${PYQT4_ROOT_DIR}/PyQt4
##
## #------ PyQt_src ------
## # nothing to do
## ##
#------ QScintilla ------
export QSCINTILLA_ROOT_DIR=${INST_ROOT}/QScintilla-2.7
exportp LD_LIBRARY_PATH ${QSCINTILLA_ROOT_DIR}/lib
##
## #------ QScintilla_src ------
## # nothing to do
## ##
#------ boost ------ 
export BOOST_ROOT_DIR=${INST_ROOT}/boost-1.52.0
exportp LD_LIBRARY_PATH ${BOOST_ROOT_DIR}/lib
##
## #------ boost_src ------
## # nothing to do
## ##
#------ Swig ------ 
export SWIG_ROOT_DIR=${INST_ROOT}/swig-2.0.8
export SWIG_LIB=${SWIG_ROOT_DIR}/share/swig/2.0.8
exportp PATH ${SWIG_ROOT_DIR}/bin
##
## #------ Swig_src ------
## # nothing to do
## ##
#------ freetype ------
export FREETYPE_ROOT_DIR=${INST_ROOT}/freetype-2.4.11
exportp LD_LIBRARY_PATH ${FREETYPE_ROOT_DIR}/lib
##
## #------ freetype_src ------
## # nothing to do
## ##
#------ freeimage ------
export FREEIMAGE_ROOT_DIR=${INST_ROOT}/freeimage-3.15.4
exportp PATH ${FREEIMAGE_ROOT_DIR}/bin 
exportp LD_LIBRARY_PATH ${FREEIMAGE_ROOT_DIR}/lib
##
## #------ freeimage_src ------
## # nothing to do
## ##
#------ cmake ------
export CMAKE_ROOT_DIR=${INST_ROOT}/cmake-2.8.10.2
exportp PATH ${CMAKE_ROOT_DIR}/bin
##
## #------ cmake_src ------
## # nothing to do
## ##
#------ gl2ps ------
export GL2PS_ROOT_DIR=${INST_ROOT}/gl2ps-1.3.8
exportp PATH ${GL2PS_ROOT_DIR}/bin 
exportp LD_LIBRARY_PATH ${GL2PS_ROOT_DIR}/lib
##
## #------ gl2ps_src ------
## # nothing to do
## ##
#------ tbb ------
export TBB_ROOT_DIR=${INST_ROOT}/tbb-30_018oss
exportp PATH ${TBB_ROOT_DIR}/bin/intel64/cc4.1.0_libc2.4_kernel2.6.16.21 
exportp LD_LIBRARY_PATH ${TBB_ROOT_DIR}/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21
##
## #------ tbb_src ------
## # nothing to do
## ##
#------ OpenCascade ------
export CAS_ROOT_DIR=${INST_ROOT}/OCCT-6.7.0
exportp PATH ${CAS_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${CAS_ROOT_DIR}/lib:${CAS_ROOT_DIR}/lin/lib
# Variable for Foundation Classes : 
export CSF_UnitsLexicon=${CAS_ROOT_DIR}/src/UnitsAPI/Lexi_Expr.dat 
export CSF_UnitsDefinition=${CAS_ROOT_DIR}/src/UnitsAPI/Units.dat 
# Variable for DataExchange : 
export CSF_SHMessage=${CAS_ROOT_DIR}/src/SHMessage
export CSF_XSMessage=${CAS_ROOT_DIR}/src/XSMessage 
# Variable for Font : 
export CSF_MDTVFontDirectory=${CAS_ROOT_DIR}/src/FontMFT 
export CSF_MDTVTexturesDirectory=${CAS_ROOT_DIR}/src/Textures 
# Activation of OCCT Kernel multithreading :
export MMGT_REENTRANT=1
# this variable only needed for DRAWEXE
export CASROOT=${CAS_ROOT_DIR}
##
## #------ OpenCascade_src ------
## # nothing to do
## ##
#------ Qwt ------
export QWT_ROOT_DIR=${INST_ROOT}/qwt-5.2.1
exportp LD_LIBRARY_PATH ${QWT_ROOT_DIR}/lib 
##
## #------ Qwt_src ------
## # nothing to do
## ##
#------ OmniORB ------
export OMNIORB_ROOT_DIR=${INST_ROOT}/omniORB-4.1.6
export OMNIORBPY_ROOT_DIR=${OMNIORB_ROOT_DIR}
exportp PYTHONPATH ${OMNIORB_ROOT_DIR}/lib:${OMNIORB_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages
exportp PATH ${OMNIORB_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${OMNIORB_ROOT_DIR}/lib 
##
## #------ OmniORB_src ------
## # nothing to do
## ##
#------ Hdf5 ------
export HDF5_ROOT_DIR=${INST_ROOT}/hdf5-1.8.10
exportp PATH ${HDF5_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${HDF5_ROOT_DIR}/lib
##
## #------ Hdf5_src ------
## # nothing to do
## ##
#------ cgnslib ------
export CGNS_ROOT_DIR=${INST_ROOT}/cgnslib-3.1.3
exportp PATH ${CGNS_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${CGNS_ROOT_DIR}/lib
##
## #------ cgnslib_src ------
## # nothing to do
## ##
#------ Med ------
export MEDFILE_ROOT_DIR=${INST_ROOT}/med-3.0.7
exportp PATH ${MEDFILE_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${MEDFILE_ROOT_DIR}/lib
##
## #------ Med_src ------
## # nothing to do
## ##
#------ Metis ------ 
export METIS_ROOT_DIR=${INST_ROOT}/metis-4.0
##
## #------ Metis_src ------
## # nothing to do
## ##
#------ Scotch ------ 
export SCOTCH_ROOT_DIR=${INST_ROOT}/scotch-5.1.11
##
## #------ Scotch_src ------
## # nothing to do
## ##
#------ ParaView ------
export PARAVIEW_ROOT_DIR=${INST_ROOT}/ParaView-3.98.1
export VTK_DIR=${PARAVIEW_ROOT_DIR}/lib/cmake/paraview-3.98
export PV_PLUGIN_PATH=${PARAVIEW_ROOT_DIR}/lib/paraview-3.98
exportp PATH ${PARAVIEW_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${PARAVIEW_ROOT_DIR}/lib/paraview-3.98
exportp PYTHONPATH ${PARAVIEW_ROOT_DIR}/lib/paraview-3.98:${PARAVIEW_ROOT_DIR}/lib/paraview-3.98/site-packages:${PARAVIEW_ROOT_DIR}/lib/paraview-3.98/site-packages/paraview
##
## #------ ParaView_src ------
## # nothing to do
## ##
#------ numpy ------
export NUMPY_ROOT_DIR=${INST_ROOT}/numpy-1.7.1
exportp PYTHONPATH ${NUMPY_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages
exportp PATH ${NUMPY_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${NUMPY_ROOT_DIR}/lib
##
## #------ numpy_src ------
## # nothing to do
## ##
#------ libBatch ------ 
export LIBBATCH_ROOT_DIR=${INST_ROOT}/libBatch-2.1.0
exportp LD_LIBRARY_PATH ${LIBBATCH_ROOT_DIR}/lib
exportp PYTHONPATH ${LIBBATCH_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages
##
## #------ libBatch_src ------
## # nothing to do
## ##
#------ expat ------
export EXPAT_ROOT_DIR=${INST_ROOT}/expat-2.0.1
exportp PATH ${EXPAT_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${EXPAT_ROOT_DIR}/lib
##
## #------ expat_src ------
## # nothing to do
## ##
#------ Graphviz ------
export GRAPHVIZ_ROOT_DIR=${INST_ROOT}/graphviz-2.30.0
exportp PATH ${GRAPHVIZ_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${GRAPHVIZ_ROOT_DIR}/lib:${GRAPHVIZ_ROOT_DIR}/lib/graphviz
##
## #------ Graphviz_src ------
## # nothing to do
## ##
#------ Doxygen ------
export DOXYGEN_ROOT_DIR=${INST_ROOT}/doxygen-1.8.3.1
exportp PATH ${DOXYGEN_ROOT_DIR}/bin
##
## #------ Doxygen_src ------
## # nothing to do
## ##
#------ Sphinx ------
export SPHINX_ROOT_DIR=${INST_ROOT}/Sphinx-1.1.3
exportp PYTHONPATH ${SPHINX_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages
exportp PATH ${SPHINX_ROOT_DIR}/bin
##
## #------ Sphinx_src ------
## # nothing to do
## ##
## #------ netgen ------
## export NETGEN_ROOT_DIR=${INST_ROOT}/netgen-4.9.13
## exportp LD_LIBRARY_PATH ${NETGEN_ROOT_DIR}/lib
## ##
## #------ netgen_src ------
## # nothing to do
## ##
## #------ MeshGems ------
## export MESHGEMS_ROOT_DIR=${INST_ROOT}/MeshGems-1.1/Products
## exportp LD_LIBRARY_PATH ${MESHGEMS_ROOT_DIR}/lib/Linux_64 
## exportp PATH ${MESHGEMS_ROOT_DIR}/bin/Linux_64 
## # license activation
## export LICENSE_FILE=/data/tmpsalome/salome/prerequis/install/LICENSE/dlim8.var.sh
## if [ -e "${LICENSE_FILE}" ] ; then
##     export SIMULOGD_LICENSE_FILE=29029@soleil
##     export DISTENE_LICENSE_FILE='Use global envvar: DLIM8VAR'
##     export DISTENE_LICENCE_FILE_FOR_MGCLEANER=${LICENSE_FILE}
##     export DISTENE_LICENCE_FILE_FOR_YAMS=${LICENSE_FILE}
##     source ${LICENSE_FILE}
## fi
## ##
## #------ MeshGems_src ------
## # nothing to do
## ##
## #------ homard ------
## export HOMARDHOME=${INST_ROOT}/homard-10.7
## exportp PATH ${HOMARDHOME}/HOMARD_V10.7_64
## ##
## #------ homard_src ------
## # nothing to do
## ##
#------ libxml2 ------ 
export LIBXML2_ROOT_DIR=${INST_ROOT}/libxml2-2.9.0
exportp PATH ${LIBXML2_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${LIBXML2_ROOT_DIR}/lib
##
## #------ libxml2_src ------
## # nothing to do
## ##
#------ wso2 ------ 
export WSO2_ROOT_DIR=${INST_ROOT}/wso2-wsf-cpp-2.1.0
exportp PATH ${WSO2_ROOT_DIR}/bin
exportp LD_LIBRARY_PATH ${WSO2_ROOT_DIR}/lib
##
## #------ wso2_src ------
## # nothing to do
## ##
#------ simanio ------ 
export SIMANIO_ROOT_DIR=${INST_ROOT}/simanio-1.0.0
exportp LD_LIBRARY_PATH ${SIMANIO_ROOT_DIR}/lib
##
## #------ simanio_src ------
## # nothing to do
## ##
#------ KERNEL ------
export KERNEL_ROOT_DIR=${INST_ROOT}/KERNEL_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${KERNEL_ROOT_DIR}/bin/salome
    exportp LD_LIBRARY_PATH ${KERNEL_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${KERNEL_ROOT_DIR}/bin/salome:${KERNEL_ROOT_DIR}/lib/salome:${KERNEL_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
  fi
fi
##
#------ KERNEL_src ------
export KERNEL_SRC_DIR=${INST_ROOT}/KERNEL_SRC_7.3.0
##
#------ GUI ------
export GUI_ROOT_DIR=${INST_ROOT}/GUI_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${GUI_ROOT_DIR}/bin/salome
    exportp LD_LIBRARY_PATH ${GUI_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${GUI_ROOT_DIR}/bin/salome:${GUI_ROOT_DIR}/lib/salome:${GUI_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
    export VTK_AUTOLOAD_PATH=${GUI_ROOT_DIR}/lib/paraview
  fi
fi
##
#------ GUI_src ------
export GUI_SRC_DIR=${INST_ROOT}/GUI_SRC_7.3.0
##
#------ GEOM ------
export GEOM_ROOT_DIR=${INST_ROOT}/GEOM_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${GEOM_ROOT_DIR}/bin/salome
    exportp LD_LIBRARY_PATH ${GEOM_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${GEOM_ROOT_DIR}/bin/salome:${GEOM_ROOT_DIR}/lib/salome:${GEOM_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
  fi
fi
##
#------ GEOM_src ------
export GEOM_SRC_DIR=${INST_ROOT}/GEOM_SRC_7.3.0
##
#------ MED ------
export MED_ROOT_DIR=${INST_ROOT}/MED_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${MED_ROOT_DIR}/bin/salome
    exportp LD_LIBRARY_PATH ${MED_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${MED_ROOT_DIR}/bin/salome:${MED_ROOT_DIR}/lib/salome:${MED_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
    export AM2CMAKE_FORCE_GENERATION=1
  fi
fi
##
#------ MED_src ------
export MED_SRC_DIR=${INST_ROOT}/MED_SRC_7.3.0
##
#------ SMESH ------
export SMESH_ROOT_DIR=${INST_ROOT}/SMESH_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${SMESH_ROOT_DIR}/bin/salome
    exportp LD_LIBRARY_PATH ${SMESH_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${SMESH_ROOT_DIR}/bin/salome:${SMESH_ROOT_DIR}/lib/salome:${SMESH_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
  fi
fi
##
#------ SMESH_src ------
export SMESH_SRC_DIR=${INST_ROOT}/SMESH_SRC_7.3.0
##
#------ PARAVIS ------
export PARAVIS_ROOT_DIR=${INST_ROOT}/PARAVIS_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${PARAVIS_ROOT_DIR}/bin/salome
    exportp LD_LIBRARY_PATH ${PARAVIS_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${PARAVIS_ROOT_DIR}/bin/salome:${PARAVIS_ROOT_DIR}/lib/salome:${PARAVIS_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
    export PV_PLUGIN_PATH="${PARAVIS_ROOT_DIR}/lib/paraview;${PV_PLUGIN_PATH}"
  fi
fi
##
#------ PARAVIS_src ------
export PARAVIS_SRC_DIR=${INST_ROOT}/PARAVIS_SRC_7.3.0
export ACCEPT_PARAVIS_WARNINGS=1
##
#------ HEXABLOCK ------
export HEXABLOCK_ROOT_DIR=${INST_ROOT}/HEXABLOCK_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${HEXABLOCK_ROOT_DIR}/bin/salome
    exportp LD_LIBRARY_PATH ${HEXABLOCK_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${HEXABLOCK_ROOT_DIR}/bin/salome:${HEXABLOCK_ROOT_DIR}/lib/salome:${HEXABLOCK_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
  fi
fi
##
#------ HEXABLOCK_src ------
export HEXABLOCK_SRC_DIR=${INST_ROOT}/HEXABLOCK_SRC_7.3.0
##
#------ HEXABLOCKPLUGIN ------
export HEXABLOCKPLUGIN_ROOT_DIR=${INST_ROOT}/HEXABLOCKPLUGIN_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp LD_LIBRARY_PATH ${HEXABLOCKPLUGIN_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${HEXABLOCKPLUGIN_ROOT_DIR}/lib/salome:${HEXABLOCKPLUGIN_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
  fi
fi
##
#------ HEXABLOCKPLUGIN_src ------
export HEXABLOCKPLUGIN_SRC_DIR=${INST_ROOT}/HEXABLOCKPLUGIN_SRC_7.3.0
##
## #------ NETGENPLUGIN ------
## export NETGENPLUGIN_ROOT_DIR=${INST_ROOT}/NETGENPLUGIN_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp LD_LIBRARY_PATH ${NETGENPLUGIN_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${NETGENPLUGIN_ROOT_DIR}/lib/salome:${NETGENPLUGIN_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ NETGENPLUGIN_src ------
## export NETGENPLUGIN_SRC_DIR=${INST_ROOT}/NETGENPLUGIN_SRC_7.3.0
## ##
## #------ GHS3DPLUGIN ------
## export GHS3DPLUGIN_ROOT_DIR=${INST_ROOT}/GHS3DPLUGIN_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp LD_LIBRARY_PATH ${GHS3DPLUGIN_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${GHS3DPLUGIN_ROOT_DIR}/lib/salome:${GHS3DPLUGIN_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ GHS3DPLUGIN_src ------
## export GHS3DPLUGIN_SRC_DIR=${INST_ROOT}/GHS3DPLUGIN_SRC_7.3.0
## ##
## #------ GHS3DPRLPLUGIN ------
## export GHS3DPRLPLUGIN_ROOT_DIR=${INST_ROOT}/GHS3DPRLPLUGIN_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${GHS3DPRLPLUGIN_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${GHS3DPRLPLUGIN_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${GHS3DPRLPLUGIN_ROOT_DIR}/bin/salome:${GHS3DPRLPLUGIN_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ GHS3DPRLPLUGIN_src ------
## export GHS3DPRLPLUGIN_SRC_DIR=${INST_ROOT}/GHS3DPRLPLUGIN_SRC_7.3.0
## ##
## #------ BLSURFPLUGIN ------
## export BLSURFPLUGIN_ROOT_DIR=${INST_ROOT}/BLSURFPLUGIN_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${BLSURFPLUGIN_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${BLSURFPLUGIN_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${BLSURFPLUGIN_ROOT_DIR}/bin/salome:${BLSURFPLUGIN_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ BLSURFPLUGIN_src ------
## export BLSURFPLUGIN_SRC_DIR=${INST_ROOT}/BLSURFPLUGIN_SRC_7.3.0
## ##
## #------ HexoticPLUGIN ------
## export HexoticPLUGIN_ROOT_DIR=${INST_ROOT}/HexoticPLUGIN_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${HexoticPLUGIN_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${HexoticPLUGIN_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${HexoticPLUGIN_ROOT_DIR}/bin/salome:${HexoticPLUGIN_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ HexoticPLUGIN_src ------
## export HexoticPLUGIN_SRC_DIR=${INST_ROOT}/HexoticPLUGIN_SRC_7.3.0
## ##
## #------ COMPONENT ------
## export COMPONENT_ROOT_DIR=${INST_ROOT}/COMPONENT_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${COMPONENT_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${COMPONENT_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${COMPONENT_ROOT_DIR}/bin/salome:${COMPONENT_ROOT_DIR}/lib/salome:${COMPONENT_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ COMPONENT_src ------
## export COMPONENT_SRC_DIR=${INST_ROOT}/COMPONENT_SRC_7.3.0
## ##
## #------ PYCALCULATOR ------
## export PYCALCULATOR_ROOT_DIR=${INST_ROOT}/PYCALCULATOR_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${PYCALCULATOR_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${PYCALCULATOR_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${PYCALCULATOR_ROOT_DIR}/bin/salome:${PYCALCULATOR_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ PYCALCULATOR_src ------
## export PYCALCULATOR_SRC_DIR=${INST_ROOT}/PYCALCULATOR_SRC_7.3.0
## ##
## #------ CALCULATOR ------
## export CALCULATOR_ROOT_DIR=${INST_ROOT}/CALCULATOR_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${CALCULATOR_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${CALCULATOR_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${CALCULATOR_ROOT_DIR}/bin/salome:${CALCULATOR_ROOT_DIR}/lib/salome:${CALCULATOR_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ CALCULATOR_src ------
## export CALCULATOR_SRC_DIR=${INST_ROOT}/CALCULATOR_SRC_7.3.0
## ##
#------ HELLO ------
export HELLO_ROOT_DIR=${INST_ROOT}/HELLO_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${HELLO_ROOT_DIR}/bin/salome
    exportp LD_LIBRARY_PATH ${HELLO_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${HELLO_ROOT_DIR}/bin/salome:${HELLO_ROOT_DIR}/lib/salome:${HELLO_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
  fi
fi
##
#------ HELLO_src ------
export HELLO_SRC_DIR=${INST_ROOT}/HELLO_SRC_7.3.0
##
#------ PYHELLO ------
export PYHELLO_ROOT_DIR=${INST_ROOT}/PYHELLO_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${PYHELLO_ROOT_DIR}/bin/salome
    exportp PYTHONPATH ${PYHELLO_ROOT_DIR}/bin/salome:${PYHELLO_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
  fi
fi
##
#------ PYHELLO_src ------
export PYHELLO_SRC_DIR=${INST_ROOT}/PYHELLO_SRC_7.3.0
##
## #------ ATOMGEN ------
## export ATOMGEN_ROOT_DIR=${INST_ROOT}/ATOMGEN_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${ATOMGEN_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${ATOMGEN_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${ATOMGEN_ROOT_DIR}/bin/salome:${ATOMGEN_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ ATOMGEN_src ------
## export ATOMGEN_SRC_DIR=${INST_ROOT}/ATOMGEN_SRC_7.3.0
## ##
## #------ ATOMSOLV ------
## export ATOMSOLV_ROOT_DIR=${INST_ROOT}/ATOMSOLV_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${ATOMSOLV_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${ATOMSOLV_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${ATOMSOLV_ROOT_DIR}/bin/salome:${ATOMSOLV_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ ATOMSOLV_src ------
## export ATOMSOLV_SRC_DIR=${INST_ROOT}/ATOMSOLV_SRC_7.3.0
## ##
## #------ ATOMIC ------
## export ATOMIC_ROOT_DIR=${INST_ROOT}/ATOMIC_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${ATOMIC_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${ATOMIC_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${ATOMIC_ROOT_DIR}/bin/salome:${ATOMIC_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ ATOMIC_src ------
## export ATOMIC_SRC_DIR=${INST_ROOT}/ATOMIC_SRC_7.3.0
## ##
#------ LIGHT ------
export LIGHT_ROOT_DIR=${INST_ROOT}/LIGHT_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp LD_LIBRARY_PATH ${LIGHT_ROOT_DIR}/lib/salome
  fi
fi
##
#------ LIGHT_src ------
export LIGHT_SRC_DIR=${INST_ROOT}/LIGHT_SRC_7.3.0
##
## #------ PYLIGHT ------
## export PYLIGHT_ROOT_DIR=${INST_ROOT}/PYLIGHT_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${PYLIGHT_ROOT_DIR}/bin/salome
##     exportp PYTHONPATH ${PYLIGHT_ROOT_DIR}/bin/salome
##   fi
## fi
## ##
## #------ PYLIGHT_src ------
## export PYLIGHT_SRC_DIR=${INST_ROOT}/PYLIGHT_SRC_7.3.0
## ##
## #------ RANDOMIZER ------
## export RANDOMIZER_ROOT_DIR=${INST_ROOT}/RANDOMIZER_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${RANDOMIZER_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${RANDOMIZER_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${RANDOMIZER_ROOT_DIR}/bin/salome:${RANDOMIZER_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ RANDOMIZER_src ------
## export RANDOMIZER_SRC_DIR=${INST_ROOT}/RANDOMIZER_SRC_7.3.0
## ##
## #------ SIERPINSKY ------
## export SIERPINSKY_ROOT_DIR=${INST_ROOT}/SIERPINSKY_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${SIERPINSKY_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${SIERPINSKY_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${SIERPINSKY_ROOT_DIR}/bin/salome:${SIERPINSKY_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ SIERPINSKY_src ------
## export SIERPINSKY_SRC_DIR=${INST_ROOT}/SIERPINSKY_SRC_7.3.0
## ##
#------ YACS ------
export YACS_ROOT_DIR=${INST_ROOT}/YACS_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PATH ${YACS_ROOT_DIR}/bin/salome
    exportp LD_LIBRARY_PATH ${YACS_ROOT_DIR}/lib/salome
    exportp PYTHONPATH ${YACS_ROOT_DIR}/bin/salome:${YACS_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
  fi
fi
##
#------ YACS_src ------
export YACS_SRC_DIR=${INST_ROOT}/YACS_SRC_7.3.0
##
#------ YACSGEN ------
export YACSGEN_ROOT_DIR=${INST_ROOT}/YACSGEN_7.3.0
if [ -n "${ENV_FOR_LAUNCH}" ] ; then
  if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
    exportp PYTHONPATH ${YACSGEN_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages
  fi
fi
##
#------ YACSGEN_src ------
export YACSGEN_SRC_DIR=${INST_ROOT}/YACSGEN_SRC_7.3.0
##
## #------ JOBMANAGER ------
## export JOBMANAGER_ROOT_DIR=${INST_ROOT}/JOBMANAGER_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${JOBMANAGER_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${JOBMANAGER_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${JOBMANAGER_ROOT_DIR}/bin/salome:${JOBMANAGER_ROOT_DIR}/lib/salome:${JOBMANAGER_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ JOBMANAGER_src ------
## export JOBMANAGER_SRC_DIR=${INST_ROOT}/JOBMANAGER_SRC_7.3.0
## ##
## #------ SAMPLES_src ------
## export DATA_DIR=${INST_ROOT}/SAMPLES_SRC_7.3.0
## ##
## #------ TUTORIAL_src ------
## ##
## #------ HOMARD ------
## export HOMARD_ROOT_DIR=${INST_ROOT}/HOMARD_7.3.0
## if [ -n "${ENV_FOR_LAUNCH}" ] ; then
##   if [ "${ENV_FOR_LAUNCH}" = "1" ] ; then
##     exportp PATH ${HOMARD_ROOT_DIR}/bin/salome
##     exportp LD_LIBRARY_PATH ${HOMARD_ROOT_DIR}/lib/salome
##     exportp PYTHONPATH ${HOMARD_ROOT_DIR}/lib/salome:${HOMARD_ROOT_DIR}/lib/python${PYTHON_VERSION}/site-packages/salome
##   fi
## fi
## ##
## #------ HOMARD_src ------
## export HOMARD_SRC_DIR=${INST_ROOT}/HOMARD_SRC_7.3.0
## ##
## #------ xdata ------ 
## export XDATAROOT=${INST_ROOT}/xdata-0.9.9
## exportp PATH ${XDATAROOT}/bin
## exportp PYTHONPATH ${XDATAROOT}/lib/python${PYTHON_VERSION}/site-packages/xdata
## ##
## #------ xdata_src ------
## # nothing to do
## ##
#------ HXX2SALOME ------
export HXX2SALOME_ROOT_DIR=${INST_ROOT}/HXX2SALOME_7.3.0/bin
exportp PATH ${HXX2SALOME_ROOT_DIR}
##
#------ HXX2SALOME_src ------
# nothing to do
##
#------ HXX2SALOMEDOC ------
#nothing to do
##
#------ HXX2SALOMEDOC_src ------
# nothing to do
##
## #------ DOCUMENTATION ------
## export DOCUMENTATION_ROOT_DIR=${INST_ROOT}/DOCUMENTATION_7.3.0
## ##
## #------ DOCUMENTATION_src ------
## export DOCUMENTATION_SRC_DIR=${INST_ROOT}/DOCUMENTATION_SRC_7.3.0
## ##
## #------ SIMAN ------
## export SIMAN_ROOT_DIR=${INST_ROOT}/SIMAN_1.0
## ##
## #------ SIMAN_src ------
## export SIMAN_SRC_DIR=${INST_ROOT}/SIMAN_SRC_1.0
## ##
