#!/bin/bash
action=$(zenity --height=350 --list --text "Choose action" --radiolist --column Pick --column Action TRUE run FALSE compile FALSE clean  FALSE "run_gencase" False "preprocess"  False "compile_gencase" FALSE clean_gencase FALSE clean_src FALSE clean_all  FALSE done);
case "$action" in
run) mpiexec -np $numproc "$appname-$runoption" | zenity --text-info --width 830   --height=950 
     source ../../action.sh
;;
compile) make | zenity --text-info --width 830   --height=950 #| (zenity --progress --title "Building FEMuS " --pulsate)
         source ../../action.sh
;;
done) 
;;
run_gencase) numproc=$(zenity --scale --text "how many proc" --min-value=1 --max-value=16 --value=2 --step 1);
mpiexec -np $numproc ../gencase/"gencase-$runoption" |zenity --text-info --width 830   --height=950
source ../../action.sh
;;
compile_gencase) cd ../gencase
make | zenity --text-info --width 830   --height=950
cd ../"$appname"
source ../../action.sh
;;
clean) make clean
source ../../action.sh
;;
clean_gencase) cd ../gencase
make clean
cd ../"$appname"
source ../../action.sh
;;
clean_src) make src_clean
source ../../action.sh
;;
clean_all) make src_clean
make resu_clean
source ../../action.sh
;;
preprocess)
source ../../preprocess.sh
;;
esac

