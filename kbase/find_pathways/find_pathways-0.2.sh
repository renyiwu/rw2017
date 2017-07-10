#!/usr/bin/bash
#Find pathways with gene in common.
# Renyi Wu, 6-30-2017
inputname="dss8wk.txt"
outputname="dss8wk-matched2.txt"
pathwayfile="pathways-dss8wk.txt"
echo "Finding pathways for file: $inputname" >"$outputname" 
echo "Based on individual pathway collection" >>"$outputname"
for gene in $(cat "$inputname" | tr "\r\n" " ") # | sed "s/\n/ /g")
do
echo "Found $gene in the following pathway(s):" >>"$outputname" 
((grep -i "$gene" "$pathwayfile" || echo "None")|cut -c1-40)  >>"$outputname"  #cut option limits output length to 40 characters
done
