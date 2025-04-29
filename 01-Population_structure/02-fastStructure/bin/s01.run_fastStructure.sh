#module load miniconda
#conda init
#conda activate faststructure #if created, else conda env create -f environment.yml

Project_name="All_data_wo_orientalis"
OUT_FOLDER="./out/"
BED_FILE_NAME="/scratch/mc10739/01-SNPs_calling_filtering/SNPs_filtering/out/s08_bis/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered"

mkdir -p ./bin/$Project_name/
mkdir -p $OUT_FOLDER$Project_name/s00.sbatch
mkdir -p $OUT_FOLDER$Project_name/s01.faststructure
mkdir -p $OUT_FOLDER$Project_name/s02.clumpak
mkdir -p $OUT_FOLDER$Project_name/s02.5.clumpak_2
mkdir -p $OUT_FOLDER$Project_name/s03.pophelper
mkdir -p $OUT_FOLDER$Project_name/s04.Plot
mkdir -p $OUT_FOLDER$Project_name/s05.membership_coefficient


python3 ./bin/mk_faststructure_jobs.py\
    -b $BED_FILE_NAME \
    -minK 2 \
    -maxK 10 \
    -r 20 \
    -o ./bin/$Project_name/faststructure_job/ \
    -f $OUT_FOLDER$Project_name/s01.faststructure

# create a list of job to submit them in an array in jubail
bash ./bin/$Project_name/faststructure_job/00.group_submits.sh

# start parallelized job
slurm_parallel_ja_submit.sh ./bin/$Project_name/faststructure_job/faststructure_job_list.txt