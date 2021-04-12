# input files
dir = "inputdirectoryname/"
indfile_main = dir + "indfilename.ind"  # IND file
indfile_admix = dir + "famfilename.fam"  # FAM file
Qfile = dir + "Qfilename.K.Q"  # Q file
Pfile = dir + "Pfilename.K.P"  # P file
Gfile = dir + "genotypefilename.raw"  # G file

# output
out_dir = "outputdirectoryname/"
outfile = out_dir + "outputfilename.txt"

# user-supplied vars
target_pops = ["Target1","Target2",...,"TargetN"]
source_pops_var  = ["VaryingSource1","VaryingSource2",...,"VaryingSourceM"]
source_pops_const = ["ConstantSource1","ConstantSource2",...,"ConstantSourceL"]
num_reps = 1000  # number of repetitions desired for bootstrap in standard error estimation
num_reps_pval = 10000  # number of repetitions desired for bootstrap in empirical p value calculations
qval = 0.01  # threshold of values of ADMIXTURE ancestral populations to be considered in the bootstrap 