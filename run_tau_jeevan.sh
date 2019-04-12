#!/bin/sh
#to setup tau in jeevan's account

source ~/.bashrc

#remove all the profile files
rm profile*

#setting up tau
export TAU_OPTIONS='-optVerbose -optRevert'
export TAU_MAKEFILE=`pwd`/../mar29/Mar29/tau/lib/Makefile.tau-mpi-pdt

#compile with tau
../mar29/Mar29/tau/bin/tau_cc.sh src.c

#set the path for tau
export PATH=`pwd`/parmetis/bin:`pwd`/../mar29/Mar29/tau/bin:$PATH


#number of processes
np=$1

time mpiexec -n $np -f hostfile ./a.out ./big_graph_dataset 13233 6474
