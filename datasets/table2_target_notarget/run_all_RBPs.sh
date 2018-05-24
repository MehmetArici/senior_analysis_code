#! /bin/bash


filename="../Table0_list_of_RBPs/SpongeWebsite_RBP_list.txt"

while read -r line
do
   RBP="$line"
   echo "python identify_target_background_generic_run_db.py ${RBP} "
   #python identify_target_background_generic_run_db.py ${RBP} 
done <  "$filename"
