#!/bin/bash
cd /data/backup/Yuni/CRC_Liver/

my_gpu_array=(2 3 4 5)
echo 'my gpu numbers:'
for(( i=0;i<${#my_gpu_array[@]};i++)) do
#${#array[@]}获取数组长度用于循环
echo ${my_gpu_array[i]};
done;

echo 'tiles batches:'
for i in {0..3}; do
t1=$[i*10];
t2=$[t1+10];

if (("$t2" >= 178))
then
	t2=177
fi

echo "with gpu:'${my_gpu_array[$i]}'";
echo "'$t1','$t2'";
{
python classification2yuni.py '${my_gpu_array[$i]}' t1 t2 
}&
done
wait
