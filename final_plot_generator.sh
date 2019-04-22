echo  "-------------------------------------------" 
echo "running AllToAll [Baseline].." 

make clean
cp fix_8.c src.c
make

rm -rf boxplot* 1>&2 2>/dev/null
rm -rf avgPlot* 1>&2 2>/dev/null

for np in `seq 8 8 64`
do
	avg=0
	for i in `seq 1 5`
	do
		echo "running_without_opt_np : "$np", iteration : "$i 
		./run_tau_jeevan.sh $np datasets/wiki-talk_eds_5021410_nds_2394385 5021410 2394385
		temp=$(pprof | tail -20 | head -1 | tr  -s " " | cut -d " " -f  3 | tr -d "," )
		echo "size : "$np" , total_time : "$temp >> boxplot_without_opt
		#avg=$(( avg+temp ))
	done
	#avg=$(( avg / 5))
	#echo "size : "$np" avg_time : "$avg >> avgPlot_without_opt
done

echo  "-------------------------------------------" 

make clean
cp fix_10_mpi_in_place.c src.c
make

echo "running with Optimizatoin .." 


for np in `seq 8 8 64`
do
	avg=0
	for i in `seq 1 5`
	do
		echo "running_with_opt_np : "$np", iteration : "$i 
		./run_tau_jeevan.sh $np  datasets/wiki-talk_eds_5021410_nds_2394385 5021410 2394385
		temp=$(pprof | tail -20 | head -1 | tr  -s " " | cut -d " " -f  3 | tr -d "," )
		echo "size : "$np" , total_time : "$temp >> boxplot_with_opt
		#avg=$(( avg+temp ))
	done
	#avg=$(( avg / 5))
	#echo "size : "$np" avg_time : "$avg >> avgPlot_with_opt
done


