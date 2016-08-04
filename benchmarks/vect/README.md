Vectorization Experiment on the Xeon Phi
========================================

This directory contains a vectorization experiment for the optimized Tersoff
potential. Running the whole experiment might take some time.

Steps to run the benchmark:

1. Create the binaries by invoking

    python generate_makefiles.py <path-to-lammps-root>

2. Run the benchmark wherever you want by calling

    ./measure-host.sh

3. Fetch results from output files via

    python enter_database.py

4. Congratulations! The data is now available either in data.dat as a usual
   data file that you can read using Excel/Python/R, or in data.db, which
   you can query as a SQLite database.

The potential files used for these benchmarks are taken from the LAMMPS
program (lammps.sandia.gov).
