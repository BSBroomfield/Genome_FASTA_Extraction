#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 1-0
#SBATCH -o /dev/null
#SBATCH -a 1-3

##Requirements 'GeneList.txt' - Column of Transcript Names (e.g., DVRF1_AAEL008492-RA) with expected length in column 2

GENE=`awk "NR==$SLURM_ARRAY_TASK_ID" GeneList.txt | awk '{print $1}'`

echo $GENE
mkdir ${GENE}_INDELS

grep -v ">" ${GENE}*.fasta | grep -E "[actg\*]+" | sort | uniq >> ${GENE}_INDELS.txt

L=`wc -l ${GENE}_INDELS.txt | awk '{print $1}'`
l=1

while [ "$l" -le "$L" ]; do

Allele=`awk "NR==$l" ${GENE}_INDELS.txt`

grep -B 1 $Allele ${GENE}*.fasta >> ${GENE}_INDELS/${GENE}_Indels.fasta
grep -B 1 "*" $GENE*.fasta >> ${GENE}_INDELS/${GENE}_Indels.fasta

l=$(( $l + 1 ))
done

rm ${GENE}_INDELS.txt
grep ">" ${GENE}_INDELS/${GENE}_Indels.fasta | cut -d ">" -f 2 | sort | uniq > ${GENE}_IndelList.txt

mkdir FASTA_NoIndels

~/Apps/seqkit grep -w 0 -n -v -f ${GENE}_IndelList.txt ${GENE}*.fasta -o FASTA_NoIndels/${GENE}_NoIndels.fasta
~/Apps/seqkit grep -w 0 -n -f ${GENE}_IndelList.txt ${GENE}*.fasta -o ${GENE}_INDELS/${GENE}_Indels.fasta

mv ${GENE}_IndelList.txt ${GENE}_INDELS/

I=`wc -l FASTA_NoIndels/${GENE}_NoIndels.fasta | awk '{print $1}'`
i=1

while [ "$i" -le "$I" ]; do

ID=`awk "NR==$i" FASTA_NoIndels/${GENE}_NoIndels.fasta`
SEQ=`awk "NR==$i+1" FASTA_NoIndels/${GENE}_NoIndels.fasta`

LEN=`echo -n $SEQ | wc -m`
EXPECTED=`awk "NR==$SLURM_ARRAY_TASK_ID" GeneList.txt | awk '{print $2}'`

if [[ "$LEN" == "$EXPECTED" ]]
then
	echo $ID >> ${GENE}_NoIndels_tmp.fasta
	echo $SEQ >> ${GENE}_NoIndels_tmp.fasta
else
	echo $ID >> ${GENE}_INDELS/${GENE}_Indels.fasta
	echo $SEQ >> ${GENE}_INDELS/${GENE}_Indels.fasta
fi

i=$(( $i + 2 ))
done

mv ${GENE}_NoIndels_tmp.fasta FASTA_NoIndels/${GENE}_NoIndels.fasta
