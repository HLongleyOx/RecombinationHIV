#!/bin/bash
#SBATCH -p short
#SBATCH --array 1-100
#SBATCH -o /well/fraser/users/hqh588/pacbio_updated/simulations_slim/logs/%j.out
#SBATCH -e /well/fraser/users/hqh588/pacbio_updated/simulations_slim/logs/%j.e


WINFILE=$1
WIN=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $WINFILE)

module load Python/3.7.4-GCCcore-8.3.0
source /well/fraser/users/hqh588/python/hivGeneral-skylake/bin/activate

for ((i = 50000; i <= 50600; i += 50)); do
    python generate_fastas.py $i $WIN $2
done



