#!/usr/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=24
#SBATCH --error=metatranscripts.%J.err
#SBATCH --output=metatranscripts.%J.out

#Example input:
#sbatch bacterial_metatranscriptomics.sh -d /media/bigdrive1/sam/Sponge_RNA_alignments/Scaffolds/HMA_only/FL2014_3 -s FL2014_3_scaffolds.fasta -b FL2014_3 -t FL2014_3.taxonomy.tsv -u FL2014_3_DNA_R1.fastq -v FL2014_3_DNA_R2.fastq -x FL2014_3_RNA_R1.fastq -y FL2014_3_RNA_R2.fastq

Help()
{
   # Display Help
   echo "Required flags"
   echo
   echo "Syntax: scriptTemplate [-d|h|s|a|b|c|e]"
   echo "-h     Display this help menu"
   echo "Required flags:"
   echo "-d     Input directory with all files"
   echo "-s     Assembly scaffolds"
   echo "-b     Base name"
   echo "-t     Contig taxonomy file"
   echo "-u     DNA forward reads"
   echo "-v     DNA reverse reads"
   echo "-x     RNA forward reads"
   echo "-y     RNA reverse reads"
   echo
}

#Define flags
while getopts "d:s:b:t:u:v:x:y:h" flag
do
    case "${flag}" in
        d) input_dir=${OPTARG};;
        s) scaffolds=${OPTARG};;
        b) basename=${OPTARG};;
        t) contig_taxonomy=${OPTARG};;
        u) dna_fwd_read=${OPTARG};;
        v) dna_rev_read=${OPTARG};;
        x) rna_fwd_read=${OPTARG};;
        y) rna_rev_read=${OPTARG};;
                h) # display Help
                        Help
                        exit;;
    esac
done

#move to directory where all files are
cd $input_dir

echo "Files in input directory"
echo
ls

echo
echo "File names as defined by user"
echo "Input directory: " $input_dir;
echo "Scaffolds: " $scaffolds;
echo "DNA reads: " $dna_fwd_read $dna_rev_read;
echo "RNA reads: " $rna_fwd_read $rna_rev_read;
echo "Contig taxonomy file: " $contig_taxonomy;
echo "All files to be prefixed with: " $basename
echo


#Confirm files defined by flags exist
if [ -f "$scaffolds" ]; then
        echo "$scaffolds found."
else 
        echo "$scaffolds does not exist."
        exit
fi

if [ -f "$dna_fwd_read" ]; then
        echo "$dna_fwd_read found."
else 
        echo "$dna_fwd_read does not exist."
        exit
fi

if [ -f "$dna_rev_read" ]; then
        echo "$dna_rev_read found."
else 
        echo "$dna_rev_read does not exist."
        exit
fi

if [ -f "$rna_fwd_read" ]; then
        echo "$rna_fwd_read found."
else 
        echo "$rna_fwd_read does not exist."
        exit
fi

if [ -f "$rna_rev_read" ]; then
        echo "$rna_rev_read found."
else 
        echo "$rna_rev_read does not exist."
        exit
fi

if [ -f "$contig_taxonomy" ]; then
        echo "$contig_taxonomy found."
else 
        echo "$contig_taxonomy does not exist."
        exit
fi

#Extract only bacterial contigs
python /home/sam/Slurm_scripts/Sponge_Project_scripts/pull_bacterial_contigs.py $input_dir $scaffolds $contig_taxonomy

echo "Bacterial contigs extracted"

#Run prodigal on scaffold file to get gff and faa files
prodigal -i ${scaffolds/.fasta/_bacterial_only.fasta} -a ${scaffolds/.fasta/_bacterial_only.faa} -f gff -p meta -o ${scaffolds/.fasta/_bacterial_only.gff}

echo "Prodigal run complete"

#Convert gff to gtf
grep -v "#" ${scaffolds/.fasta/_bacterial_only.gff}| grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,4,5,7,9 |  awk -v OFS='\t' '{print $1,"PROKKA","CDS",$2,$3,".",$4,".","gene_id " $5}' \
> ${scaffolds/.fasta/_bacterial_only.gtf}

echo "Prokka gff file converted to gtf file"

#Get mapped reads and transcripts per gene
bbmap.sh \
in1=$dna_fwd_read \
in2=$dna_rev_read \
trimclip \
outm=${scaffolds/.fasta/_DNA.bam} \
ref=$scaffolds

echo "DNA reads mapped to all contigs"

samtools sort -n -O bam -o ${scaffolds/.fasta/_sorted_DNA.bam} ${scaffolds/.fasta/_DNA.bam}

echo "Mapped DNA reads sorted"

featureCounts -T 10 -t CDS -g gene_id -F GTF -p -B -P -d 50 -D 600 -a ${scaffolds/.fasta/_bacterial_only.gtf} -o ${scaffolds/.fasta/_DNA_count.txt} ${scaffolds/.fasta/_sorted_DNA.bam}

echo "DNA reads counted per bacterial gene"

bbmap.sh \
in1=$rna_fwd_read \
in2=$rna_rev_read \
trimclip \
outm=${scaffolds/.fasta/_RNA.bam} \
ref=$scaffolds

echo "RNA reads mapped to all contigs"

samtools sort -n -O bam -o ${scaffolds/.fasta/_sorted_RNA.bam} ${scaffolds/.fasta/_RNA.bam}

echo "Mapped RNA reads sorted"

featureCounts -T 10 -t CDS -g gene_id -F GTF -p -B -P -d 50 -D 600 -a ${scaffolds/.fasta/_bacterial_only.gtf} -o ${scaffolds/.fasta/_RNA_count.txt} ${scaffolds/.fasta/_sorted_RNA.bam}

echo "RNA transcripts counted per bacterial gene"

rm *.bam
rm *.summary

#Get rid of headers and replace spaces with single tab in count files
sed -e '1,2d' < ${scaffolds/.fasta/_RNA_count.txt} > ${scaffolds/.fasta/_RNA_count.tab}
rm ${scaffolds/.fasta/_RNA_count.txt} 
awk -v OFS='\t' '{ $1=$1; print }' ${scaffolds/.fasta/_RNA_count.tab} > ${scaffolds/.fasta/_RNA_count.txt}
rm ${scaffolds/.fasta/_RNA_count.tab}

sed -e '1,2d' < ${scaffolds/.fasta/_DNA_count.txt} > ${scaffolds/.fasta/_DNA_count.tab}
rm ${scaffolds/.fasta/_DNA_count.txt}
awk -v OFS='\t' '{ $1=$1; print }' ${scaffolds/.fasta/_DNA_count.tab} > ${scaffolds/.fasta/_DNA_count.txt}
rm ${scaffolds/.fasta/_DNA_count.tab}

echo "Count tables formatted"

#Run kofamscan on all bacterial genes
eval "$(conda shell.bash hook)"
conda activate kofamscan
/home/sam/Tools/kofamscan/kofamscan ${scaffolds/.fasta/_bacterial_only.faa} --cpu 16 -o ${scaffolds/.fasta/_bacterial_only_kegg_output_detailed} -p /home/sam/Tools/kofamscan/db/profiles -k /home/sam/Tools/kofamscan/db/ko_list -f detail-tsv
conda deactivate 

echo "KofamScan run complete"

#Extract only reliable KEGG results
python /home/sam/Slurm_scripts/Sponge_Project_scripts/strict_kegg_extractor.py $input_dir

echo "Reliable KOfamScan annotations identified and stored"

#Get COG annots
eval "$(conda shell.bash hook)"
conda activate eggnog
emapper.py -i ${scaffolds/.fasta/_bacterial_only.faa} -o ${scaffolds/.fasta/_eggNOG_output} --output_dir $input_dir --cpu 16 --data_dir /media/bigdrive1/Databases/eggnog_mapper_latest
conda deactivate 

#Do TPM calculations
python /home/sam/Useful_scripts/bacterial_final_file_generator.py $input_dir $basename
