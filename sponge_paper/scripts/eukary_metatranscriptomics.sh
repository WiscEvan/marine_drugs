#!/usr/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=24
#SBATCH --error=metatranscripts.%J.err
#SBATCH --output=metatranscripts.%J.out

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
   echo "-g     gff file from AUGUSTUS"
   echo "-f     fa file from AUGUSTUS"   
   echo "-t     Contig taxonomy file"
   echo "-u     DNA forward reads"
   echo "-v     DNA reverse reads"
   echo "-x     RNA forward reads"
   echo "-y     RNA reverse reads"
   echo
}

#Define flags
while getopts "d:s:b:g:f:t:u:v:x:y:h" flag
do
    case "${flag}" in
        d) input_dir=${OPTARG};;
        s) scaffolds=${OPTARG};;
        b) basename=${OPTARG};;
        g) gff_file=${OPTARG};;
        f) fa_file=${OPTARG};;
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
echo "AUGUSTUS gff file: " $gff_file;
echo "AUGUSTUS fa file: " $fa_file;
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


#Convert gff to gtf
grep -v "#" $gff_file | cut -f1 -d ';' | grep "gene"| cut -f1,4,5,7,9| awk -v OFS='\t' '{print $1,"AUGUSTUS","CDS",$2,$3,".",$4,".","gene_id " $5}' > ${gff_file/.gff/.gtf}

#echo "AUGUSTUS gff file converted to gtf file"

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

featureCounts -T 10 -t CDS -g gene_id -F GTF -p -B -P -d 50 -D 600 -a ${gff_file/.gff/.gtf} -o ${scaffolds/.fasta/_DNA_count.txt} ${scaffolds/.fasta/_sorted_DNA.bam}

echo "DNA reads counted per eukaryotic gene"

bbmap.sh \
in1=$rna_fwd_read \
in2=$rna_rev_read \
trimclip \
outm=${scaffolds/.fasta/_RNA.bam} \
ref=$scaffolds

echo "RNA reads mapped to all contigs"

samtools sort -n -O bam -o ${scaffolds/.fasta/_sorted_RNA.bam} ${scaffolds/.fasta/_RNA.bam}

echo "Mapped RNA reads sorted"

featureCounts -T 10 -t CDS -g gene_id -F GTF -p -B -P -d 50 -D 600 -a ${gff_file/.gff/.gtf} -o ${scaffolds/.fasta/_RNA_count.txt} ${scaffolds/.fasta/_sorted_RNA.bam}

echo "RNA transcripts counted per eukaryotic gene"

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

echo "Count tables formatted correctly"

#Run kofamscan on all sponge genes
eval "$(conda shell.bash hook)"
conda activate kofamscan
/home/sam/Tools/kofamscan/kofamscan $fa_file --cpu 16 -o ${fa_file/.fa/kegg_output_detailed} -p /home/sam/Tools/kofamscan/db/profiles -k /home/sam/Tools/kofamscan/db/ko_list -f detail
conda deactivate 

echo "KofamScan run complete"

#Extract only reliable KEGG results
python /home/sam/Slurm_scripts/Sponge_Project_scripts/eukary_strict_kegg_extractor.py $input_dir

echo "Reliable KOfamScan annotations identified and stored"

#Get COG annots
eval "$(conda shell.bash hook)"
conda activate eggnog
emapper.py -i $fa_file -o ${fa_file/.fa/_eggNOG_output} --output_dir $input_dir --cpu 16 --data_dir /media/bigdrive1/Databases/eggnog_mapper_latest
conda deactivate 

echo "eggNOG annotations identified and stored"

#Do TPM calculations
python /home/sam/Useful_scripts/eukary_file_generator.py $input_dir $basename
