# input files
#dir = "inputdirectoryname/"
dir = "example/"
#indfile_main = dir + "indfilename.ind"  # IND file
indfile_main = dir + "example_ind_file.ind"  # IND file
#indfile_admix = dir + "famfilename.fam"  # FAM file
indfile_admix = dir + "example_file.fam"  # FAM file
#Qfile = dir + "Qfilename.K.Q"  # Q file
Qfile = dir + "example_file.6.Q"  # Q file
#Pfile = dir + "Pfilename.K.P"  # P file
Pfile = dir + "example_file.6.P"  # P file
#Gfile = dir + "genotypefilename.raw"  # G file
Gfile = dir + "example_geno_file.raw"  # G file

# output
#out_dir = "outputdirectoryname/"
out_dir = "example/"
#outfile = out_dir + "outputfilename.txt"
outfile = out_dir + "example_results.txt"

# user-supplied vars
#target_pops = ["Target1","Target2",...,"TargetN"]
target_pops = ["Som50Eng50"]
#source_pops_var  = ["VaryingSource1","VaryingSource2",...,"VaryingSourceM"]
source_pops_var  = ["Somali", "Spanish", "Iranian"]
#source_pops_const = ["ConstantSource1","ConstantSource2",...,"ConstantSourceL"]
source_pops_const = ["English"]
num_reps = 1000  # number of repetitions desired for bootstrap in standard error estimation
num_reps_pval = 10000  # number of repetitions desired for bootstrap in empirical p value calculations
qval = 0.01  # threshold of values of ADMIXTURE ancestral populations to be considered in the bootstrap
