#! /usr/bin/env python
import os
from module_generator import Generator,Module,Service,HXX2SALOMEComponent
#
# grab from environment variables of interest
kernel_root_dir=os.environ["KERNEL_ROOT_DIR"]
gui_root_dir=os.environ["GUI_ROOT_DIR"]
yacs_root_dir=os.environ["YACS_ROOT_DIR"]
salome_env="/homesd/msandro/software/salome/env_products.sh"

# create a context, that provide environment information 
context={'update':1,
         "makeflags":"-j2",
         "prerequisites":salome_env,
         "kernel":kernel_root_dir,
         "gui":gui_root_dir,
         "yacs":yacs_root_dir
        }
#
# specify a SALOME modules called MATHEMATICS, including one SALOME component called MATH
module="MATHEMATICS"
components_name=["MATH"]
components=[]
for compo_name in components_name:
  # cpp root dir for the cpp library= /dsk4/dsk4/salome/BUILD/MYAPPLICATION/MATHCPP
  cpp_root_dir="/homesd/msandro/software/femus/contrib/Myapplication/MATHCPP"
  # print
  print compo_name,   components_name,    cpp_root_dir
  # generation of class components
  components.append(HXX2SALOMEComponent("MATH.hxx", "libMATHCXX.so",cpp_root_dir ) )

# we ask to install the module in /dsk4/dsk4/salome/MATHEMATICS
module_root_dir="/homesd/msandro/software/salome/MATHEMATICS"

# generate and compile the salome component with yacsgen
g=Generator(Module(module,components=components,prefix=module_root_dir),context)
g.generate()
g.bootstrap()
g.configure()
g.make()
g.install()

# generate a SALOME application containing :
#   - the MATHEMATICS module, 
#   - the native modules GUI, YACS
g.make_appli("appli_math", restrict=["KERNEL"], 
                      altmodules={"GUI":yacs_root_dir, 
                                  "YACS":yacs_root_dir})
