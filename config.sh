module load openmpi/1.5.4
module load numpy/1.6.1
module load scipy/0.10.0
module load biopython/1.61
module load EMBOSS/6.4.0

export WGAC=/net/eichler/vol2/eee_shared/assemblies/hg19/wgac/genomicSuperDup.tab
export SUNKS=/net/eichler/vol2/eee_shared/assemblies/hg19/sunks/hg19_sunks.bed.gz
export REFERENCE_DIR=/net/eichler/vol2/eee_shared/assemblies/hg19
export REFERENCE=${REFERENCE_DIR}/ucsc.hg19.fasta
export REPEATS=/net/eichler/vol2/eee_shared/assemblies/hg19/repeats/repeats_and_trf.bed
