#!bin/bash
 mpiexec  -np 8 femlcore_1-opt
 mkdir b20
 mv resu b20
 mv DATA/parameters.in b20
 mv DATA/parameters.in2 DATA/parameters.in
 echo programma completato| mailx -A gmail -s "fine_corsa" danycerr@gmail.com
 echo programma completato| mailx -A gmail -s "fine_corsa" simone.bna@gmail.com
 mpiexec  -np 8 femlcore_1-opt
 mkdir b60
 mv resu b60
 mv DATA/parameters.in b60
 echo programma completato 2 | mailx -A gmail -s "fine_corsa2" danycerr@gmail.com
 echo programma completato 2| mailx -A gmail -s "fine_corsa2" simone.bna@gmail.com