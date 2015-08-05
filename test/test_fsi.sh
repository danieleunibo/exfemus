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

APP="fsi_1"
APP_NAME=${APP/%_*/""}
APP_VER=${APP/#*_/""}
echo  "APP_NAME" $APP_NAME "APP_VER" $APP_VER
export FM_MYAPP=$APP
export METHOD=opt
export PETSC_ARCH=linux-opt 
cd ../
# ======================================================================
# ==== GO TO USER_APPL AND COPY OR CREATE THE APPLICATION FM_MYAPP =====
# ======================================================================
cp -r template_appl/$APP_NAME/   USER_APPL/$APP
cd USER_APPL/$FM_MYAPP
echo "APPLICATION=" $FM_MYAPP "APP_NAME=" $APP_NAME "APP_VER=" $APP_VER "METHOD=" $METHOD "in USER_APPL/"$FM_MYAPP
echo  "path " $PWD 
mkdir RESU
make gencase
mpiexec -np 1 ../gencase/gencase-opt
make
mpiexec -np 1 $FM_MYAPP-opt