
# CIRM project on isolated points #

Code associated to the paper "Towards a classification of isolated j-invariants"
by Abbey Bourdon, Sachi Hashimoto, Timo Keller, Zev Klagsbrun, David Lowry-Duda, Travis Morrison, Filip Najman, Himanshu Shukla

System requirements: Magma version at least 2.27 is required to read the data files.

Input:

On a specific curve:     j-invariant
On a data set of curves: txt file with list of j-invariants
                         delimited by newline character, e.g. adelicgenusgt0curves.txt

Output:

On a specific curve:     [j-invariant, empty list] or
                         [j-invariant, set of pairs (level and degree)]

The list [< n_1,d_1>, ..., < n_k,d_k>] contains all points such that any
isolated point $x \in X_1(N)$ with $j(x)= j$ maps down under the natural
projection map to an isolated point of degree $d_i$ on $X_1(n_i)$ for
$1 \leq i \leq k$.


## INSTRUCTIONS ##

1. This repository depends David Roe's modified version of David Zywina's
   repository found here <https://github.com/roed314/OpenImage.git> as well as <https://github.com/AndrewVSutherland/ell-adic-galois-images> by by Jeremy Rouse, Andrew V. Sutherland, and David Zureick-Brown. To clone these repos, run

   ```
   git clone https://github.com/roed314/OpenImage.git
   git clone https://github.com/AndrewVSutherland/ell-adic-galois-images.git
   ```

   To run for a specific curve given by a j-invariant j, do the following.

2. Open magma.

3. Attach the required files for the above repositories, by running

   ```
   AttachSpec("path/to/OpenImage/OpenImage.spec");
   Attach("path/to/ell-adic-galois-images/groups/gl2.m");
   ```

4. Attach the spec file for this repository by running

   ```
   AttachSpec("/path/to/isolated.spec");
   ```

5. Define

   ```
   path := OpenImageContext("path/to/data-files")
   ```

   where the path leads to the data files in the OpenImage repository.

6. Run the command

   ```
   NotIsolated(j, path);
   ```


## Advanced instructions ##

To run for the entire list of curves in the file adelicgenusgt0curves.txt you also need to
install GNU parallel. Then execute the following command

    parallel --joblog JOB_LOG_FILE_NAME --shuf --timeout TIMEOUT \
      --eta -jNUMBER_OF_CORES \
      magma -b seq:={} doline.m ::: {LINE_NUMBERS_OF_INTEREST} >> OUTPUT_FILE_NAME.out

Details of the last command:

- TIMEOUT: Set it to the amount of time in seconds for which each instance of
  NotIsolated.m has to run e.g. 200

- NUMBER_OF_CORES: Set it to the number of cores which have to engaged during
  parallel execution. Check the number of cores available or the machine
  beforehand.

- JOB_LOG_FILE_NAME: Set it to the desired file name for job logs. An exitval
  value of -1 for a particular instance in this file will indicate if a
  particular instance was terminated because it required more time than
  TIMEOUT.

- LINE_NUMBERS_OF_INTEREST: Set it to the line numbers corresponding to the
  curves in the file Curves.txt for which NotSporadic.m has to be executed.
  E.g. 11, 13, 19,... or one can use a continuous range as in magma like
  1..100.

- OUTPUT_FILE_NAME: Set it to the desired file name for the output file. This
  file contains output of every successful execution i.e. all the executions
  which were not timed out.


Note: Users with a Bash terminal can run install_dependencies.sh using the
command ./install_dependencies.sh. This will install GNU parallel 2023 version
and also the OpenImage github repo.


## License Information ##

Copyright (c) 2024, held by the isolated points team: Bourdon, Hashimoto, Keller, Klagsbrun, Lowry-Duda, Morrison, Najman, and Shukla.

The magma files

- appendix.m
- Code_for_X1(17)_and_X1(24).m
- doline.m
- extrafilters.m
- isolatedpoints.m
- testlozanorobledo.m

are made available as documentation, code, and educational materials under the
Creative Commons Attribution 4.0 (CC-BY-4.0) license. See
https://creativecommons.org/licenses/by/4.0/ for more details, or
[LICENSE](./LICENSE) for a copy of the license text.

We provide these in the hope that they will be useful. Please contact us if
there are licensing reasons preventing you from making use of this code. We note
that [this license is compatible with the GPL](https://www.gnu.org/licenses/license-list.en.html#ccby)
as long as you make proper attribution and follow the other guidelines of CC-BY-4.0.


The paper *Towards a classification of isolated j-invariants* (including early
versions of the tex files) is not covered by this license.
