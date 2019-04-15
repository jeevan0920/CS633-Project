echo  "-------------------------------------------" >> output_Isend_10times
echo "running AllToAll [Baseline].." >> ouptut_Isend_10times
cp fix_8.c src.c

avg_Isend=0
for i in `seq 1 10`
do
	./run_tau_siva.sh 10 datasets/wiki-talk_eds_5021410_nds_2394385 5021410 2394385
	temp=$(pprof | tail -20 | head -1 | tr  -s " " | cut -d " " -f  3 | tr -d "," )
	#echo "temp = "$temp	
	avg_Isend=$(( $avg_Isend + $temp ))
	#echo "avg =  "$avg_Isend
done

avg_Isend=`expr $avg_Isend / 10`
echo  "average allToall Baseline: "$avg_Isend >> output_Isend_10times


echo  "-------------------------------------------" >> output_Isend_10times
cp fix_10_mpi_in_place.c src.c
echo "running with Optimizatoin .." >> output_Isend_10times
avg_Isend=0
for i in `seq 1 10`
do
	./run_tau_siva.sh 10 datasets/wiki-talk_eds_5021410_nds_2394385 5021410 2394385
	temp=$(pprof | tail -20 | head -1 | tr  -s " " | cut -d " " -f  3 | tr -d "," )
	#echo "temp = "$temp	
	avg_Isend=$(( $avg_Isend + $temp ))
	#echo "avg =  "$avg_Isend
done

avg_Isend=`expr $avg_Isend / 10`
echo  "average Isend and deallocation and mpi_in_place: "$avg_Isend >> output_Isend_10times


