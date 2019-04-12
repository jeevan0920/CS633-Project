#!/bin/sh
#to setup tau in jeevan's account

source ~/.bashrc

#setting up tau
export TAU_OPTIONS='-optVerbose -optRevert'
export TAU_MAKEFILE=`pwd`/../mar29/Mar29/tau/lib/Makefile.tau-mpi-pdt
../mar29/Mar29/tau/bin/tau_cc.sh
export PATH=`pwd`/parmetis/bin:`pwd`/../mar29/Mar29/tau/bin:$PATH


time mpiexec -n 10 -f hostfile ./a.out ./datasets/big_graph_slash_dot 516575 77357
