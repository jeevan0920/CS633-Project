echo  "-------------------------------------------" 
echo "running AllToAll [Baseline].." 

cp fix_8.c src.c

rm -rf columnPlot* 1>&2 2>/dev/null

for np in `seq 8 8 64`
do
		echo "running_without_opt_np : "$np 
		./run_tau_jeevan.sh $np datasets/wiki-talk_eds_5021410_nds_2394385 5021410 2394385
		temp=$(pprof | tail -20 | head -1 | tr  -s " " | cut -d " " -f  3 | tr -d "," )
		allgatherv=$(pprof -s | grep "MPI_Allgatherv()"  | head -1  | tr -s " "  | cut -d " " -f 2)
		alltoallv=$(pprof -s | grep "MPI_Alltoallv()"  | head -1  | tr -s " "  | cut -d " " -f 2)
		echo "size : "$np" , allgatherv : "$allgatherv" , alltoallv : "$alltoallv >> columnPlot_without_opt
done

echo  "-------------------------------------------" 

cp fix_10_mpi_in_place.c src.c

echo "running with Optimizatoin .." 


for np in `seq 8 8 64`
do
		echo "running_with_opt_np : "$np 
		./run_tau_jeevan.sh $np datasets/wiki-talk_eds_5021410_nds_2394385 5021410 2394385
		temp=$(pprof | tail -20 | head -1 | tr  -s " " | cut -d " " -f  3 | tr -d "," )
		allgatherv=$(pprof -s | grep "MPI_Allgatherv()"  | head -1  | tr -s " "  | cut -d " " -f 2)
		mpiisend=$(pprof -s | grep "MPI_Isend()"  | head -1  | tr -s " "  | cut -d " " -f 2)
		mpiirecv=$(pprof -s | grep "MPI_Irecv()"  | head -1  | tr -s " "  | cut -d " " -f 2)
		sum=$(( mpiisend + mpiirecv ))
		echo "size : "$np" , allgatherv : "$allgatherv" , MPI_Isend() + MPI_Recv() : "$sum >> columnPlot_with_opt
done




