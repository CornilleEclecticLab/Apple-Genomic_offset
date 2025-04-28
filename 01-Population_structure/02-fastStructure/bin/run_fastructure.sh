#module load miniconda
#conda init
#conda activate faststructure #if created, else conda env create -f environment.yml
mkdir -p ./out/s00.sbatch
mkdir -p ./out/s01.faststructure

python3 ./bin/mk_faststructure_jobs.py\
    -b /scratch/mc10739/01-SNPs_calling_filtering/SNPs_filtering/out/gatk/s08.pop_gen_filter/chr_concat.geno20_maf005.anno.syno.LDpruned \
    -minK 2 \
    -maxK 10 \
    -r 20 \
    -o ./bin/faststructure_job/

# create a list of job to submit them in an array in jubail
bash ./bin/bin/faststructure_job/00launch_submits.sh

# start parallelized job
slurm_parallel_ja_submit.sh bin/faststructure_job_list.txt