Compilation Instructions
========================
	 # goto code directory
	 > cd code 

	 run the make file
	 > make


Executing source code
========================
	 > mpiexec -f hostfile -np <number-processes> ./src.x <data-set-path>  <number-edges> <number-vertices>


Profiling with tau
=======================
	# set environment variables for tau compiler
	> export TAU_OPTIONS='-optVerbose -optRevert'
	> export TAU_MAKEFILE=`pwd`/tau/lib/Makefile.tau-mpi-pdt
	
	# set the path for tau compiler
	> export PATH=<path-to-tau-bin>:$PATH
	
	# compiling with tau
	>tau/bin/tau_cc.sh src.c 

Executing tau-compiled code
===========================
	# execute likely normally, on success it will generate profiling files
	> mpiexec -np <number-processes> -f hostfile ./a.out <data-set-path>  <number-edges> <number-vertices>
	
	# read the profiling files using pprof or paraprof
	> pprof  
	> paraprof


