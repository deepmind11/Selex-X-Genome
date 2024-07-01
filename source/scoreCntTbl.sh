#!/bin/bash


count_table=$1
fit_id=$2
base_dr=$(dirname $count_table)
filename_with_ext=$(basename $count_table)
filename="${filename_with_ext%.*}"
echo $filename

#Creating seq.txt file from the count table.
awk '{print $1}' $count_table > $base_dr/$filename.txt

#Score the seq.txt file using ProboundTools
java -cp /burg/home/hg2604/bin/ProBoundTools/ProBoundTools/target/ProBound-jar-with-dependencies.jar  proBoundTools/App -c "loadMotifCentralModel($fit_id).buildConsensusModel().addNScoring().inputTXT($base_dr/$filename.txt).bindingModeScores(/dev/stdout)" | awk '{print $2}' >  ${base_dr}/${filename}_${fit_id}.txt

#Deleting the temp files
rm $base_dr/$filename.txt
