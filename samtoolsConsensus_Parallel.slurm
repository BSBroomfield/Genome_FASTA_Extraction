#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 0-1
#SBATCH -o %x.o%j
#SBATCH -a 1

module load apps/gcc/samtools/1.21
module load apps/bioinf
module load apps/seqtk/1.0/gcc-4.8.5

#TO RUN N GENES IN PARALLEL, SET '#SBATCH -a 1-X' WHERE X=number of genes in 'GeneList.txt'

#REQUIREMENTS IN DIRECTORY
#"SampleList.txt" : File containing Sample ID prefix for BAM files, one sample ID per line
#"GeneList.txt" : Two column file containing GFF Gene ID in column 1 (case sensitive) and a gene name in column 2 (can be custom names)

##SET PATHWAYS
GFF=~/scratch/WholeGenome_data/Ref_Genome/VectorBase-68_AaegyptiLVP_AGWG.gff #Path to GFF
BAMPATH=~/scratch/WholeGenome_data/Indel_Realigned_Bam #Path to BAM directory (Note: Omit final '/')

##SET PARAMETERS
MINCOV=98	#Minimum mean coverage of a gene to include a sample
MINDP=6		#Minimum mean depth of a gene to include a sample and minimum site depth to include sites in consensus sequence
MINMQ=30	#Minimum mapping quality to include sites in consensus
MINBQ=20	#Minimum base quality to include sites in consensus
N_LIM=0.001	#Maximum proportion of N-bases allowed in final sequence (e.g., 0.05 = final sequnce has <=5% N bases)
INDELS=yes	#Option to toggle indels on or off. 'yes' includes indels in final sequence 'no' ignores indels

n=$SLURM_ARRAY_TASK_ID

GENE_ID=`awk "NR==$n" GeneList.txt | awk '{print $1}'`
GENE_NAME=`awk "NR==$n" GeneList.txt | awk '{print $2}'`

mkdir CodingRegions
mkdir CodingRegions/${GENE_NAME}
mkdir ${GENE_NAME}

#EXTRACT UNIQUE TRANSCRIPTS FOR EACH GENE ID
grep ${GENE_ID} ${GFF} | grep mRNA | awk '{print $9}' | cut -d ";" -f 1 | cut -d "=" -f 2 > ${GENE_NAME}/${GENE_NAME}_TranscriptIDs.txt

T=`wc -l ${GENE_NAME}/${GENE_NAME}_TranscriptIDs.txt | awk '{print $1}'`
t=1
while [ "$t" -le "$T" ]; do

TRANSCRIPT=`awk "NR==$t" ${GENE_NAME}/${GENE_NAME}_TranscriptIDs.txt`
mkdir ${GENE_NAME}/${TRANSCRIPT}

grep $TRANSCRIPT $GFF | grep CDS | awk '{print $1":"$4"-"$5}' > CodingRegions/${GENE_NAME}/${TRANSCRIPT}_CDS.txt

K=`wc -l SampleList.txt | awk '{print $1}'`
k=1
while [ "$k" -le "$K" ]; do

ID=`awk "NR==$k" SampleList.txt`

##CALCULATE COVERAGE FOR EACH EXON
E=`wc -l CodingRegions/${GENE_NAME}/${TRANSCRIPT}_CDS.txt | awk '{print $1}'`
e=1
while [ "$e" -le "$E" ]; do
REGION=`awk "NR==$e" CodingRegions/${GENE_NAME}/${TRANSCRIPT}_CDS.txt`

samtools coverage -q ${MINMQ} -Q ${MINBQ} -r $REGION ${BAMPATH}/${ID}*.bam | awk 'FNR==2' >> ${GENE_NAME}/${TRANSCRIPT}/${ID}_${GENE_NAME}_COV.txt

e=$(($e + 1))
done

#AVERAGE EXON COVERAGE STATISTICS AND ONLY KEEP SAMPLES PASSING THRESHOLDS
MEANCOV=`awk -v N=6 '{ sum += $N } END { if (NR > 0) print sum / NR }'  ${GENE_NAME}/${TRANSCRIPT}/${ID}_${GENE_NAME}_COV.txt | cut -d "." -f 1`
MEANDP=`awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'  ${GENE_NAME}/${TRANSCRIPT}/${ID}_${GENE_NAME}_COV.txt | cut -d "." -f 1`
MEANBQ=`awk -v N=8 '{ sum += $N } END { if (NR > 0) print sum / NR }'  ${GENE_NAME}/${TRANSCRIPT}/${ID}_${GENE_NAME}_COV.txt | cut -d "." -f 1`
MEANMQ=`awk -v N=9 '{ sum += $N } END { if (NR > 0) print sum / NR }'  ${GENE_NAME}/${TRANSCRIPT}/${ID}_${GENE_NAME}_COV.txt | cut -d "." -f 1`

echo "${ID} ${MEANCOV} ${MEANDP} ${MEANBQ} ${MEANMQ}" >> ${GENE_NAME}/${TRANSCRIPT}/${GENE_NAME}_${TRANSCRIPT}_Coverage.txt
rm ${GENE_NAME}/${TRANSCRIPT}/${ID}_${GENE_NAME}_COV.txt

if (( "$MEANCOV" >= "$MINCOV" && "$MEANDP" >= "$MINDP" ))
then

echo "${ID}" >> ${GENE_NAME}/${TRANSCRIPT}/${GENE_NAME}_${TRANSCRIPT}_SampleList.txt

fi

k=$(($k+1))
done

sed -i '1 i\Sample_ID Mean_Cov Mean_DP Mean_BQ Mean_MQ' ${GENE_NAME}/${TRANSCRIPT}/${GENE_NAME}_${TRANSCRIPT}_Coverage.txt

S=`wc -l ${GENE_NAME}/${TRANSCRIPT}/${GENE_NAME}_${TRANSCRIPT}_SampleList.txt | awk '{print $1}'`
s=1
while [ "$s" -le "$S" ]; do
SAMPLE=`awk "NR==$s" ${GENE_NAME}/${TRANSCRIPT}/${GENE_NAME}_${TRANSCRIPT}_SampleList.txt`

##EXTRACT CONSENSUS FOR EACH EXON
E=`wc -l CodingRegions/${GENE_NAME}/${TRANSCRIPT}_CDS.txt | awk '{print $1}'`
e=1
while [ "$e" -le "$E" ]; do
REGION=`awk "NR==$e" CodingRegions/${GENE_NAME}/${TRANSCRIPT}_CDS.txt`

samtools consensus -r $REGION -f FASTA -o ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}_Exon${e}.fasta \
--min-MQ $MINMQ --min-BQ $MINBQ \
--show-del $INDELS --show-ins $INDELS \
-A -d $MINDP \
${BAMPATH}/${SAMPLE}*.bam

e=$(($e + 1))
done

##CONCATENATE EXONS INTO SINGLE FASTA FORMAT FILE
awk 'FNR>=2' ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}_Exon* | tr -d '\n' | awk -v header=${SAMPLE}_${GENE_NAME}_${TRANSCRIPT} 'FNR==1 {print ">"header}1' > ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}.fasta
rm ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}_Exon*

##REVERSE COMPLEMENT GENES ON NEGATIVE STRAND
STRAND=`grep $GENE_ID $GFF | grep CDS | awk 'FNR==1 {print $7}'`

if [[ "$STRAND" == "-" ]]
then

seqtk seq -r ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}.fasta > ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}_rev.fasta
mv ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}_rev.fasta ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}.fasta

fi

##FINAL 'N' QC
SeqLen=`awk "NR==2" ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}.fasta | wc -m`
SeqN=`awk "NR==2" ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}.fasta | sed 's/[^N]//g' -`
NumN=`echo -n $SeqN | wc -m`
PropN=`awk "BEGIN {printf \"%.4f\", $NumN / $SeqLen}"`
echo $SeqLen $SeqN $NumN $PropN

if echo "$N_LIM >= $PropN" | bc -l | grep -q 1
then

cat ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}.fasta >> ${GENE_NAME}/${TRANSCRIPT}/${GENE_NAME}_${TRANSCRIPT}.fasta
echo ${SAMPLE} >> ${GENE_NAME}/${TRANSCRIPT}/Final_${GENE_NAME}_${TRANSCRIPT}_SampleIDs.txt

fi

rm ${GENE_NAME}/${TRANSCRIPT}/${SAMPLE}_${GENE_NAME}_${TRANSCRIPT}.fasta

s=$(($s + 1))
done

t=$(($t + 1))
done
