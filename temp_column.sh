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
                echo "size : "$np" , allgatherv : "$allgatherv" , MPI_Isend() : "$mpiisend",  MPI_Recv() : "$mpiirecv >> columnPlot_with_opt
done
