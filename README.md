# LINADMIX

LINADMIX is an algorithm that tests how well a mixture of ancient populations explains the genetics of a modern population.

It consists of 4 Python scripts: run_linadmix.py, the master script; config.py, which contains the user-defined variables; admix.py, which defines the Admix class and its methods; and linadmix_tools.py, which contains helper functions.

The master script (run_linadmix.py) expects 5 input files:

-the .fam and .raw files from PLINK

-the .P and .Q files from ADMIXTURE

-a .ind file in the EIGENSTRAT format

The input file names and their directory should be specified in the config file. In the same script, the user can specify the name of the output directory and file, as well as the source and target populations, the desired number of bootstrap repetitions, and the minimum Q-value for the variance calculations.

The output is a text file that gives the target population and the source populations, followed by the mixing coefficients and standard errors (corresponding to the source populations) and the P value.

The script takes no command line inputs and can be run using the command **python run_linadmix.py**.

The *example* directory contains 7 files. Five are examples of the input files LINADMIX expects. The config file shows how the variables are set. The .txt
file is the output using the example inputs.