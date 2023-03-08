Data files required:

1) gon_data.m (Gonality data for $X_1(N)$ when $N\geq 300$)
2) Curves.txt (file containing the list of representative curves for
   j-invariants in LMFDB which do not correspond to adelic genus 0 and do not have
   CM)

Input:

On a specific curve:     Cremona label/Weierstrass coefficient array.
On a data set of curves: txt file with list of Cremona labels/Weierstrass coefficient array
                         delimited by newline character, e.g. Curves.txt

Output: 
On a specific curve:     [j-invariant, Weierstrass coefficient array, true] or 
                         [j-invariant, Weierstrass coefficient array, false, set of
                          pairs (label and degree)]
On a set of curves:      Joblog file with name JOB_LOG_FILE_NAME in the parallel
                         command and an output file with the name OUTPUT_FILE_NAME.

INSTRUCTIONS:

1)  Execute install_dependencies.sh using the command
    ./install_dependencies.sh . This will install parallel 2023 version and also
    David Roe's github repo on modular curves. This is important to run
    FindOpenImage() as an intrinsic magma function. The directory to install has
    to be entered in the format /home/....

To run for a specific curve given by E where E is a Cremona label or Weierstrass
coefficient vector execute the following commands

2)  Open magma.

3)  Execute load "NotSporadic.m";

4)  Execute NotSporadic(E);

To run for a list of curves in the file Curves.txt execute the following command

2) parallel --joblog JOB_LOG_FILE_NAME --shuf --timeout TIMEOUT --eta -jNUMBER_OF_CORES magma -b seq:={} NotSporadic.m ::: {LINE_NUMBERS_OF_INTEREST} >> OUTPUT_FILE_NAME.out

Details of the last command: 

TIMEOUT: Set it to the amount of time in
seconds for which each instance of NotSporadic.m has to run e.g. 200

NUMBER_OF_CORES: Set it to the number of cores which have to engaged during
parallel execution. Check the number of cores available or the machine
beforehand.

JOB_LOG_FILE_NAME: Set it to the desired file name for job logs. An exitval
value of -1 for a particular instance in this file will indicate if a 
particular instance was terminated because it required more time than TIMEOUT.

LINE_NUMBERS_OF_INTEREST: Set it to the line numbers corresponding to the curves
in the file Curves.txt for which NotSporadic.m has to be executed. E.g. 11, 13,
19,... or one can use a continuous range as in magma like 1..100.

OUTPUT_FILE_NAME: Set it to the desired file name for the output file. This file
contains output of every successful execution i.e. all the executions which were
not timed out.  

