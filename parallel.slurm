#!/bin/bash

# FIRST LINES TO SPECIFFY SLURM DETAILS OF THE JOB

#SBATCH --job-name=y_call
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=standard
#SBATCH --array=0-23%5
#SBATCH --output=jobs/y_call_%a.o

# Define the beginnig and the end of the range to make Y haplogroup calls (chunks of 10):

start=$((SLURM_ARRAY_TASK_ID * 10)) # Begin at 0, 10, 20, 30... 190, 200
end=$((start + 10)) # End at 10, 20, 30... 200, 210,

# Activate conda
source /mnt/primevo/shared_data/software/miniconda3/25.1.1/etc/profile.d/conda.sh 
conda activate y_call

# Add information for the duration of the execution
echo "Job $SLURM_ARRAY_TASK_ID start: $(date)"

start_time=$(date +%s) # store moment of start.

# Run the Python program
#python make_call.py --bam_list "./benchmarking_Narasimhan.tsv" --index_i $start --index_f $end --snps data/input/snps_hg38_YBrowse.csv --rg "hg38" --create_network Y --transitions Y --haplogr_tree data/input/YTree.csv --branches data/input/hap_info.csv --mismatches data/input/mm_v3.tsv

#python make_call.py --bam_list "/mnt/archgen/users/hringbauer/git/EPIDEMIC/output/tables/autoeager_seqfull_anno.tsv" --index_i $start --index_f $end --snps "data/input/YBrowse_snps_hg38-runtime.csv" --rg "hg38" --transitions "Y" --create_network Y --width 120 --height 90

python make_call.py --bam_list "/mnt/archgen/users/hringbauer/git/EPIDEMIC/output/tables/autoeager_seqfull_anno.tsv" --index_i $start --index_f $end --rg "hg38" --create_network Y --transitions Y --database YFull

end_time=$(date +%s) # store moment of end.
duration=$((end_time - start_time)) # difference between start and end as the time of execution. 

echo -e "\nJob $SLURM_ARRAY_TASK_ID end: $(date)"
echo -e "\nDuration (seconds): $duration"
