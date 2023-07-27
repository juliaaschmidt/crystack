

for file in *.res; do 
n=$(echo $file | grep -o -E '[0-9]+')
echo $n >> input_list.txt
done
