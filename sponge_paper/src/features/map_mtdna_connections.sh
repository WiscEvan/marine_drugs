#!/bin/bash

export PATH=$PATH:/home/erees/bin:/home/erees/scripts

rsync -avzh /mnt/gluster/erees/spades_input/BS-WH5_Whole*R* ./
rsync -avzh /mnt/gluster/erees/bsimplexV3.10_assembly.fna ./

defaultAssembly="bsimplexV3.10_assembly.fna"
defaultName="bsimplexV3.10ScaffoldConnections"
f_reads="BS-WH5_Whole_R1_PE_f.fastq.gz"
r_reads="BS-WH5_Whole_R2_PE_f.fastq.gz"
threads=10

#1a. Build a bowtie2 database for use in bowtie2 alignment program
./bowtie2-2.3.4.3-linux-x86_64/bowtie2-build $defaultAssembly $defaultName --threads $threads

#1b. Align the comparison dataset reads to your assembly with bowtie2.
./bowtie2-2.3.4.3-linux-x86_64/bowtie2 -x $defaultName \
  -1 $f_reads -2 $r_reads \
  -q --phred33 --very-sensitive --no-unal \
  -S $defaultName".sam" \
  --threads $threads \
  --all

#2. Convert the SAM file to a sorted BAM file , and create a BAM index.
./bin/samtools view -bS $defaultName".sam" | ./bin/samtools sort -o $defaultName".sort"
./bin/samtools view -bS $defaultName".sam" > $defaultName".bam"
#3. Tabulate the average coverage of each contig.
./scripts/fasta_length_table.pl $defaultAssembly > $defaultName".len.tab"
./bedtools2/bin/genomeCoverageBed -ibam $defaultName".sort" -g $defaultName".len.tab" > $defaultName".cov.txt"
./scripts/contig_coverage_from_bedtools.pl $defaultName".cov.txt" > $defaultName".avg.cov.txt"

#4. Calculate the average read length to be used as a parameter for the cytoscape network visualization
default_average_read_length=`./scripts/average_read_length_from_sam.py $defaultName".sam"`

#5. Run build cytoscape interaction table with attributes from
#: Albertsen et. al. Paper Nature Biotechnology 2013
#: http://madsalbertsen.github.io/multi-metagenome/docs/step10.html
bp_within_each_paired_end=500
min_scaffold_len=3000

echo "Average Read Length:"$default_average_read_length
echo "bases within each paired end cutoff:"$bp_within_each_paired_end
echo "Min Scaffold Length:"$min_scaffold_len

perl ./MadsAlbertsen-multi-metagenome-62e7eee/cytoscapeviz/cytoscapeviz.pl -i $defaultName".sam" \
  -e $bp_within_each_paired_end \
  -a $default_average_read_length \
  -m $min_scaffold_len \
	-c

tar -czvf /mnt/gluster/erees/bsimplexV3.10_connections.tar.gz $defaultName* cytoscape* condensed*
