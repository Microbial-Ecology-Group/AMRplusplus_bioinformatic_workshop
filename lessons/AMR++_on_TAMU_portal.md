
module load Nextflow

cd Desktop/

cd AMRplusplus_beta/

nextflow run main_AMR++.nf -profile singularity_workshop

cd test_results

cd multiQC

nextflow run main_AMR++.nf -profile singularity_workshop -with-singularity /home/training/qiime2.sif