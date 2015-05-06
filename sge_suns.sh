#! /bin/bash
#$ -S /bin/bash
#$ -cwd

set -e

# Send an email when the script begins, ends, aborts, or suspends.
#$ -m beas

#$ -l disk_free=20.0G
#$ -l h_vmem=6G

if [[ "$#" -ne "3" ]]
then
    echo "Usage: $0 <wgac alignments> <wgac_sequences.fasta> <output_dir>"
    exit 1
fi

INPUT=$1
REFERENCE=$2
OUTPUT_DIR=$3
REFERENCE_DIR=/tmp/wgac
REFERENCE_NAME=`basename ${REFERENCE}`

echo "Input: ${INPUT}"
echo "Reference: ${REFERENCE}"
echo "Output: ${OUTPUT_DIR}"

export P4_RSHCOMMAND=/usr/bin/rsh
export PRINT_SEQUENCES=2
export JOBDIR=$TMPDIR/work
export WORKING_DIR=$SGE_O_WORKDIR

# Load required modules.
. /etc/profile.d/modules.sh

if test ! -z $MODULESHOME; then
   module load modules modules-init/prod modules-gs/prod
   module load openmpi/1.5.4
fi

module load python/2.7.2
module load numpy/1.6.1
module load scipy/0.10.0
module load biopython/1.61
module load EMBOSS/6.4.0

# Copy reference to all nodes.
mpirun -x PATH -x LD_LIBRARY_PATH \
  --prefix $MPIBASE -mca plm ^rshd \
  -mca btl ^openib python /net/eichler/vol4/home/jlhudd/src/rsync_mpi/batch_node_copy.py \
  --source "${REFERENCE}" --dest "${REFERENCE_DIR}" \
  --pre_sync_commands "mkdir -p ${REFERENCE_DIR}" \
  --post_sync_commands="chgrp -R eichlerlab ${REFERENCE_DIR}; chmod -R g+rwx ${REFERENCE_DIR}" \
  --rsync_options "-arzv --bwlimit=20000" \
  --max_nodes 5

# Run SUN finder.
mpirun -x PATH -x LD_LIBRARY_PATH \
  --prefix $MPIBASE -mca plm ^rshd \
  -mca btl ^openib python /net/eichler/vol4/home/jlhudd/pipelines/sunks/suns/find_suns.py \
  $INPUT ${REFERENCE_DIR}/${REFERENCE_NAME} ${OUTPUT_DIR}
