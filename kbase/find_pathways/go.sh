#cat aomdss8wk.txt | while read line
#do
#grep "$line" pathways-aomdss8wk.txt
#done

for i in $(cat aomdss8wk.txt)
do
echo "$i"
grep -i "$i" pathways-aomdss8wk2.txt
#echo "\tmark"
done >log.txt
j="tnf"
for j in $(cat aomdss8wk.txt)
do
echo "$j was found in:"
grep -i $j pathways-all.txt >log.txt
done
#do
#echo $LINE
#done
