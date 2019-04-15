cp fix_9_Isend.c src.c
./jobscript.sh
echo "running Isend .." >> output_allgatherv
avg_Isend=0
for i in `seq 4 4 64`
do
	./run_tau_siva.sh $i datasets/wiki-talk_eds_5021410_nds_2394385 5021410 2394385
	temp=$( pprof | tail -19 | head -1 | tr  -s " " | cut -d " " -f  2 )
	echo "process = "$i",temp = "$temp >> output_allgatherv
done

./jobscript_8.sh
echo "running Isend . with 8 ppn." >> output_allgatherv
avg_Isend=0
for i in `seq 4 4 64`
do
        ./run_tau_siva.sh $i datasets/wiki-talk_eds_5021410_nds_2394385 5021410 2394385
        temp=$( pprof | tail -19 | head -1 | tr  -s " " | cut -d " " -f  2 )
        echo "process = "$i",ppn=8,temp = "$temp >> output_allgatherv
done

