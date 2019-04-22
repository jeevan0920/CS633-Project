for np in `seq 48 8 64`
do
        avg=0
        for i in `seq 1 5`
        do
                ./run_tau_jeevan.sh $np datasets/wiki-talk_eds_5021410_nds_2394385 5021410 2394385
                temp=$(pprof | tail -20 | head -1 | tr  -s " " | cut -d " " -f  3 | tr -d "," )
                echo "size : "$np" , total_time : "$temp >> boxplot_with_opt
                #avg=$(( avg+temp ))
        done
        #avg=$(( avg / 5))
        #echo "size : "$np" avg_time : "$avg >> avgPlot_with_opt
done
