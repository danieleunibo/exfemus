#set environment
source ../env.sh

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

APP=vofsi_1
APP_NAME=${APP/%_*/""}
APP_VER=${APP/#*_/""}
export FM_MYAPP=$APP
######   $2 ->METHOD  #################
export METHOD=opt
export PETSC_ARCH=linux-opt  
cd ../
# ======================================================================
# ==== GO TO USER_APPL AND COPY OR CREATE THE APPLICATION FM_MYAPP =====
# ======================================================================
#  create  USER_APPL/$FM_MYAPP directory
cp -r template_appl/$APP_NAME/   USER_APPL/$FM_MYAPP
cd USER_APPL/$FM_MYAPP  
echo "APPLICATION=" $FM_MYAPP "APP_NAME=" $APP_NAME "APP_VER=" $APP_VER "METHOD=" $METHOD "in USER_APPL/"$FM_MYAPP
echo  "path " $PWD
cd ../gencase
make clean
make -j4
mpiexec -np 4 gencase-opt
cd ../$FM_MYAPP
make src_clean
make -j4
mpiexec -np 4 $FM_MYAPP-opt