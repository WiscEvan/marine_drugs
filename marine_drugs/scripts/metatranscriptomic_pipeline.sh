#!/usr/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=24
#SBATCH --error=gene_summ.%J.err
#SBATCH --output=gene_summ.%J.out
#SBATCH --mail-user=samche42@gmail.com
#SBATCH --mail-type=ALL

cd /media/bigdrive1/sam/Sponge_RNA_alignments/FL2014_9

for file in `ls *.gff`; do gffread ${file} -T -o ${file/.gff/.gtf};done

for file in `ls *.faa`
do
	diamond blastp -d /media/bigdrive1/Databases/COG_data/COG_diamond.dmnd \
	-q ${file} -k 1 --max-hsps 1 \
	--outfmt 6 qseqid stitle pident evalue qlen slen \
	-o ${file/.faa/.tab}
done

python /home/sam/Useful_scripts/COG_summary.py /media/bigdrive1/sam/Sponge_RNA_alignments/FL2014_9

for bin in `ls *.fna`
do
        bbmap.sh \
        in1=FL20-9_CTTGTA_L005_polyat-clean_PE1.fastq.gz \
        in2=FL20-9_CTTGTA_L005_polyat-clean_PE2.fastq.gz  \
        outm=${bin/.fna/.bam} \
        ref=${bin}

        featureCounts -T 10 -t CDS -g gene_id -a ${bin/.fna/.gtf} -o ${bin/.fna/_count.txt} ${bin/.fna/.bam}

        rm ${bin/.fna/.bam}
done

python /home/sam/Useful_scripts/calculate_TPM_per_bin.py /media/bigdrive1/sam/Sponge_RNA_alignments/FL2014_9
