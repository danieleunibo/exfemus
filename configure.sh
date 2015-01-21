export INSTALLATION_DIR=/home/daniele/software 
export OPENMPI_NAME=openmpi-1.8.3 
export PETSC_NAME=petsc-3.5.2 
export HDF5_NAME=salome_7.4/hdf5-1.8.10 
export LIBMESH_NAME=libmesh-0.9.4 
export SALOME_NAME=salome_7.4 
export med3_NAME=med-3.0.7 
export MED_NAME=MED_7.4.0

########## Over here there MUST be the definition of       #############
########## $INSTALLATION_DIR, given when installing femus; #############
########## otherwise, write it on your own                 #############
#######################################################################


#for the GUI, check this variable
export SAL_VER=6.5.0



# ============================================================
# ================ HELP ======================================
# ============================================================

# In order to use femus run source configure.sh from femus directory.
# 
# $$ source configure.sh APPNAME_number [options]
# 
# APPNAME must be an existing application in the template_appl directory,
# else you can create an empty application to write on your own.
# Open always three tabs for your application, one in dbg mode, one 
# in opt mode and the last with the correct gencase.




# ============================================================
# =============== CHECKING APPLICATION =======================
# ============================================================

APP=$1
APP_NAME=${APP/%_*/""}
APP_VER=${APP/#*_/""}
echo  "APP_NAME" $APP_NAME "APP_VER" $APP_VER
if test "$APP_NAME" = "$1";  then APP_NAME="None";fi
echo  "APP_NAME" $APP_NAME "APP_VER" $APP_VER
cd template_appl 
ls
if [ -d  $APP_NAME ]; then
  echo $APP_NAME " is an application template - ok  "
  cd ..
else
  echo  " "$1" name is invalid: not VALID application template or _n is missing "
  cd ..
    echo
    echo " Use: source configure.sh application_n [option]"
    echo
    echo "   applications: applications_n (example nse_1, fsi_2, ..; n=version) "
    echo " 	nse (example nse_1, nse_2, ..; n=version) "
    echo " 	fsi  (example fsi_1,  fsi_2,  ..; n=version) "
    echo " 	femlcore (example femlcore_1, femlcore_2, ..; n=version) "
    echo " 	turb (example turb_1, turb_2, ..; n=version) "
    echo
    echo "   [option] : None, opt, dbg, gui, nogui, gencase "
    echo "		None = no gui in opt mode "
    echo "		opt  = salome gui in opt mode"
    echo "		dbg  = no gui in dbg mode"
    echo "		nogui= no gui in opt mode"
    echo "		gencase= no gui opt mode in gencase directory"
    echo
fi


  export FM_MYAPP=$1
######   $2 ->METHOD  #################
  OPTION=$2
  if test "$OPTION" = "dbg"; then
   export METHOD=dbg
   OPTION="nogui"
   export PETSC_ARCH=linux-dbg
  else
   export METHOD=opt
   export PETSC_ARCH=linux-opt   
  fi

  
  
# ============================================================
# =================  PATH for LIBRARIES  =====================
# ============================================================
export FEMUS_DIR=$INSTALLATION_DIR/femus/
# -----------  LIBMESH ---------------
export LIBMESH_PATH=$INSTALLATION_DIR/$LIBMESH_NAME/
# -----------  PETSC   ---------------
export PETSC_DIR=$INSTALLATION_DIR/$PETSC_NAME/
# -----------  OPENMPI  ---------------
export MPI_BIN_PATH=$INSTALLATION_DIR/$OPENMPI_NAME/bin/
export MPI_LIB_PATH=$INSTALLATION_DIR/$OPENMPI_NAME/lib64/
# -----------  HDF5  ---------------
export HDF5_PATH=$INSTALLATION_DIR/$HDF5_NAME/
# -----------  med  ---------------
export med_PATH=$INSTALLATION_DIR/$SALOME_NAME/$med3_NAME
# -----------  MED  ---------------
export MED_PATH=$INSTALLATION_DIR/$SALOME_NAME/$MED_NAME
# -----------  GUI  ---------------
export FM_GUI=gui
GUI_FEMUS="../../"$FM_GUI
GUI_CMD="code_femus"

##############   LIBRARIES ###################
######### MPI ##########
export PATH=$MPI_BIN_PATH:$PATH
export LD_LIBRARY_PATH=$MPI_LIB_PATH:$LD_LIBRARY_PATH
######## FEMUS #########
export LD_LIBRARY_PATH=$LIBMESH_PATH/lib64/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$PETSC_DIR/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HDF5_PATH/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$med_PATH/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$MED_PATH/lib/salome:$LD_LIBRARY_PATH
################################################





# ======================================================================
# ==== GO TO USER_APPL AND COPY OR CREATE THE APPLICATION FM_MYAPP =====
# ======================================================================
#  exist or  create  USER_APPL/$FM_MYAPP directory
if [ -d  USER_APPL/$FM_MYAPP ]; then
  # USER_APPL/$FM_MYAPP  exists
  echo "USER_APPL/"$FM_MYAPP " directory exists  " 
else
  # USER_APPL/$FM_MYAPP does not exist
  echo "USER_APPL/"$FM_MYAPP " directory does not  exist  "
  echo "Do you want to create the USER_APPL/"$FM_MYAPP  "application ? [yes/no]"
  read test_app
  if test "$test_app" = "yes"; then 
    cp -r template_appl/$APP_NAME/   USER_APPL/$1
  else
    return
  fi
fi
  cd USER_APPL/$FM_MYAPP
  
  echo "APPLICATION=" $FM_MYAPP "APP_NAME=" $APP_NAME "APP_VER=" $APP_VER "METHOD=" $METHOD "in USER_APPL/"$FM_MYAPP
  echo  "path " $PWD 
  
##############   OPTIONS  ###################

  if test "$OPTION" = "nogui"; then
   return;
  fi

  if test "$OPTION" = "gencase"; then
   $PWD
   cd ../gencase
   return;
  fi

  if test "$OPTION" = "gui"; then
   cd $GUI_FEMUS
   python bin/$GUI_CMD gui &
   cd ../USER_APPL/$FM_MYAPP
   echo " path " USER_APPL/$FM_MYAPP

   return;
  fi
   
   # salome style GUI ciao
   cd $GUI_FEMUS
   python bin/$GUI_CMD  salome &
   echo $PWD
   cd ../USER_APPL/$FM_MYAPP
