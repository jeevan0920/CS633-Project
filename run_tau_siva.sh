#!/bin/sh
#to setup tau in jeevan's account 

if [ $# -ne 4 ]
then 
	echo "given #arguments : "$#" , expected #arguments : 5"
	echo "Usage :: script_name np dataset #edges #vertices"
	exit
fi


source ~/.bashrc

#remove all the profile files
rm profile*

#setting up tau
export TAU_OPTIONS='-optVerbose -optRevert'
export TAU_MAKEFILE=`pwd`/../../mar29/Mar29/tau/lib/Makefile.tau-mpi-pdt

#compile with tau
../../mar29/Mar29/tau/bin/tau_cc.sh src.c

#set the path for tau
export PATH=`pwd`../../mar29/Mar29/parmetis/bin:`pwd`/../../mar29/Mar29/tau/bin:$PATH


#number of processes
np=$1
datasetPath=$2
edges=$3
verts=$4


time mpiexec -n $np -f hostfile ./a.out $datasetPath  $edges $verts
