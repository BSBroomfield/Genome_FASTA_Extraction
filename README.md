# Genome_FASTA_Extraction

Provided are scripts used in [ref] to pull FASTA sequences of candidate genes from BAM files and clean indels/truncated alleles.
All scripts were written in bash for a SLURM HPC cluster

Code is provided as is, but is open source for use and manipulation 

samtoolsConsesnus_parallel.slurm
This script does the initial sequenece extractions from the genome

Requirements:

  **Input data**
  
   Binary alignment map (BAM) files
   
   General Feature Format (GFF) file for the refernce genome used in BAM file
    
  **Software:**
  
   Samtools v1.21 or later (https://github.com/samtools/samtools)
   
   seqtk1.0 or later (https://github.com/lh3/seqtk)

  **Supplementary files** - see example files
  
   SampleList.txt = File containing Sample ID prefix for each BAM file to be analysed, one sample ID per line
   
   GeneList.txt =Two column file containing GFF Gene ID in column 1 (case sensitive) and a gene name in column 2 (can be custom names)

  **Parameter settings**
  
  MINCOV=[int]	#Minimum mean coverage of a gene to include a sample e.g., 98 = coverage >=98
  
  MINDP=[int]		#Minimum mean depth of a gene to include a sample and minimum site depth to include sites in consensus sequence e.g., 6 = minimum depth 6
  
  MINMQ=[int]	  #Minimum mapping quality to include sites in consensus e.g., 30 =  MQ>=30
  
  MINBQ=[int]  	#Minimum base quality to include sites in consensus e.g., 20 = BQ>=20
  
  N_LIM=[float]	#Maximum proportion of N-bases allowed in final sequence e.g., 0.05 = final sequnce has <=5% N bases

 **Output files**
 
A directory is created for each gene in GeneList.txt in a new directory 'CodingRegions'

Within each gene directory, a directory is made for each transcript (if a gene has multiple transcript annotations)

Each transcript has 4 output files:

1) *_Coverage.txt - gives Mean_Cov Mean_DP Mean_BQ Mean_MQ for each Sample_ID listed in SampleList.txt
3) *.fasta  - FASTA output for samples passing parameter thresholds
4) *_SampleList.txt
