#!/bin/bash
femshfile=$(zenity --file-selection  --title configure.sh);
export femshpath=`echo $femshfile | sed 's/configure.sh//'`
echo $femshpath
echo $femshfile
cd $femshpath
runoption=$(zenity --list --text "Choose method" --radiolist --column Pick --column Option TRUE opt FALSE dbg);
export appname=""
#appname=$(zenity --list --text "Choose Application name" --radiolist --column Pick --column Option TRUE main_1 FALSE femlcore FALSE fs1);
appname=$( ls USER_APPL/| zenity --list   --column appname );
if [[ -z $appname ]]
then
export appname=`zenity --entry --text "nome nuova applicazione" --entry-text ""`
fi
if [[ "$runoption" == opt ]]
then
echo source configure.sh $appname nogui
source configure.sh $appname nogui
else
source configure.sh $appname $runoption
fi
source ../../action.sh
