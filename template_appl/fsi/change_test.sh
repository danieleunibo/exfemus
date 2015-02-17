#!/bin/bash
if [ "$1" = "" ]; then
  echo "you need to select a test"
  echo "use as change_test 1"
else
  case $1 in
    0)
#         sed -i s/.*define.*DIMENSION.*(.*).*/#define\ DIMENSION\ \(2\)/ ./DATA/Domain_conf.h
# 	sed -i 's/.*define.*MATBC_INTERFACE.*/\/\/ #define\ MATBC_INTERFACE/' ./DATA/Domain_conf.h
# 	sed -i s/.*nolevels.*/\ nolevels\ \ \ \ \ \ 2/ DATA/parameters.in
# 	sed -i s/.*libmesh_gen.*/\ libmesh_gen\ \ \ 1/ DATA/parameters.in
	echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	echo ">>>               BATHE 2D fsi test                 >>>"
	echo ">>>     changed Domain_conf and parameters.in       >>>"
	echo ">>>         **    libmesh generator                 >>>"
        echo ">>>         **        2levels                       >>>"
	echo ">>>   Compile both gencase and your application     >>>"
	echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	;;
    1)
        sed -i s/.*define.*DIMENSION.*(.*).*/#define\ DIMENSION\ \(3\)/ DATA/Domain_conf.h
	sed -i 's/.*define.*MATBC_INTERFACE.*/\ #define\ MATBC_INTERFACE/' DATA/Domain_conf.h
	sed -i s/.*nolevels.*/\ nolevels\ \ \ \ \ \ 1/ DATA/parameters.in
	sed -i s/.*libmesh_gen.*/\ libmesh_gen\ \ \ 0/ DATA/parameters.in
	echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	echo ">>>               BATHE 3D fsi test                 >>>"
	echo ">>>     changed Domain_conf and parameters.in       >>>"
	echo ">>>           **    medmem mesh                     >>>"
        echo ">>>           **    graphics BC                     >>>"
	echo ">>>  now you can compile both gencase and your app  >>>"
	echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	;;
    *)
	echo "X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X"
	echo "X~X~X     NOT IMPLEMENTED YET    X~X~X"
	echo "X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X~X"
  esac
fi