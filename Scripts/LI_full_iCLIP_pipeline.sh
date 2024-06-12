
#!/bin/bash


cd /home2/2459810a/iCLIP_analysis/Hackathon_PURAB/
mkdir -p Data/Merged

cat \
/home2/2459810a/iCLIP_analysis/Hackathon_PURAB/Input_Data/Sequencing_data/iClip2_Mock_200221_S0_L001_R1_001.fastq.gz \
/home2/2459810a/iCLIP_analysis/Hackathon_PURAB/Input_Data/Sequencing_data/iClip2_Mock_200221_S0_L002_R1_001.fastq.gz \
/home2/2459810a/iCLIP_analysis/Hackathon_PURAB/Input_Data/Sequencing_data/iClip2_Mock_200221_S0_L003_R1_001.fastq.gz \
/home2/2459810a/iCLIP_analysis/Hackathon_PURAB/Input_Data/Sequencing_data/iClip2_Mock_200221_S0_L004_R1_001.fastq.gz >\
Data/Merged/PURAB_iCLIP_RAW.fastq.gz 


echo "Remove PhiX reads"

/software/bbmap-v38.90/bbduk.sh \
in=Data/Merged/PURAB_iCLIP_RAW.fastq.gz  \
out=Data/Merged/PURAB_iCLIP_rmPhiX.fastq.gz \
ref=/software/bbmap-v38.90/resources/phix174_ill.ref.fa.gz \
k=31 hdist=1 -Xmx16g threads=8 \
stats=Data/Merged/PURAB_iCLIP_rmPhiX.stats.txt

echo "Full dataset sequencing QC"

mkdir -p Data/Merged/FastQC

for file in Data/Merged/*_rmPhiX.fastq.gz
do
echo "Processing file $file"
fastqc $file -o Data/Merged/FastQC
echo ""
done 

mkdir -p Plots/QC

multiqc -f \
Data/Merged/FastQC/*_fastqc.zip \
-o Plots/QC \
-n 01_Sequencing-runs_MultiQC-report

echo "Count barcode frequencies"

zcat Data/Merged/XRN1_iCLIP_rmPhiX.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len)}' | sort | uniq -c | sort -k1,1rn>barcode_freqs

echo "Demultiplex reads"

mkdir -p Trimming/Demultiplex

java -jar /software/je_1.2/je_1.2_bundle.jar demultiplex \
F1=Data/Merged/PURAB_iCLIP_rmPhiX.fastq.gz \
BF=XRN1_barcodes_umi.txt \
RCHAR=':' \
O=Trimming/Demultiplex \
UF1=Unassigned_XRN1_iCLIP_merged_jemultiplexer.fastq \
M=XRN1_iCLIP_merged_jemultiplexer_out_stats.txt \
FASTQ_FILE_EXTENSION=fastq

echo "Demultiplexed reads QC"

mkdir -p Trimming/Demultiplex/FastQC

for file in Trimming/Demultiplex/*.fastq.gz
do
echo "Processing file $file"
fastqc $file -o Trimming/Demultiplex/FastQC
echo ""
done 

multiqc -f \
Trimming/Demultiplex/FastQC/*_fastqc.zip \
-o Plots/QC \
-n 02_Libraries_MultiQC-report

echo "Adapter trimming"

mkdir -p Trimming/Adapter

for file in Trimming/Demultiplex/*.fastq.gz
do
[[ $file == *Unassigned*.fastq.gz ]] && continue 
echo "" 
echo "Processing file $file" 
base=$(basename "$file") 
outfile="$(echo ${base} |sed -e 's/.fastq.gz/_trimmed.fastq.gz/')" 
outshort="$(echo ${base} |sed -e 's/.fastq.gz/_trimmed.fastq.tooshort.gz/')" 
outinfo="$(echo ${base} |sed -e 's/.fastq.gz/_trimmed.info/')" 
cutadapt $file \
-a AGATCGGAAGAGCGGTTCAG \
-j 4 -e 0.1 -O 1 --nextseq-trim 10 --minimum-length 15 \
-o Trimming/Adapter/$outfile \
--too-short-output Trimming/Adapter/$outshort > \
Trimming/Adapter/$outinfo 
echo "" 
done

echo "Adapter trimming QC"

multiqc -f \
Trimming/Adapter/*_trimmed.info \
-o Plots/QC \
-n 03_Cutadapt_MultiQC-report

echo "Assess rRNA contamination"

mkdir -p Trimming/rRNA

for file in Trimming/Adapter/*.fastq.gz
do
echo ""
base=$(basename $file .fastq.gz)
echo "Processing $file"
/software/bbmap-v38.90/bbduk.sh \
-Xmx24G \
in=$file \
out=Trimming/rRNA/$(echo ${base}).minus.rRNA.fastq.gz \
ref=GENOMEDIR/Hsap_rDNA_U13369.1.fa \
k=25 \
stats=Trimming/rRNA/$(echo ${base}).rRNA.stat.txt \
hdist=1
echo ""
done

echo "Unzip fastqs"


for file in Trimming/Adapter/*_trimmed.fastq.gz 
do 
gunzip -f $file 
done



#!/bin/bash

echo "Align reads"

mkdir -p Alignment 

for file in Trimming/Adapter/*_trimmed.fastq 
do 
echo "" 
echo "Processing file $file" 
base=$(basename "$file") 
outdir="$(echo ${base} |sed -e 's/_trimmed.fastq//')" 
mkdir Alignment/$outdir 
STAR --runMode alignReads --runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--genomeDir GENOMEDIR/STAR \
--sjdbGTFfile GENOMEDIR/Homo_sapiens.GRCh38.106_SINV.gtf \
--readFilesIn $file --outFileNamePrefix Alignment/$outdir/ \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 999 \
--outFilterMultimapNmax 1 \
--outFilterMismatchNoverLmax 0.04 \
--outSJfilterReads Unique \
--alignEndsType EndToEnd
echo "" 
done 


multiqc -f \
Alignment/*/*Log.final.out \
-o ./ \
-n 05_STAR_MultiQC-report


mkdir -p Alignment_idxstats

for dir in $(ls -d Alignment/*/) 
do
echo "" 
file=$dir/Aligned.sortedByCoord.out.bam 
base=$(basename "$dir") 
output="Alignment_idxstats/$(echo ${base}).idxstats.txt"  
echo "Processing file $file" 
samtools index $file 
samtools idxstats $file > $output 
done

echo "Alignment QC"

multiqc -f \
Alignment_idxstats/*idxstats.txt \
-o ./ \
-n 06_Chromosome-mapping-pre-dedup_MultiQC-report

echo "Alignment QC"

mkdir -p Dedup 

for dir in $(ls -d Alignment/*/) 
do
echo "" 
input=$dir/Aligned.sortedByCoord.out.bam 
base=$(basename "$dir") 
output="Dedup/$(echo ${base}).dedup.je.bam" 
outlog="Dedup/$(echo ${base}).dedup.je.log" 
metrics="Dedup/$(echo ${base}).dedup.je.metrics" 
java -jar /software/je_1.2/je_1.2_bundle.jar markdupes \
I=$input \
O=$output \
M=$metrics \
REMOVE_DUPLICATES=True \
MM=1 >\
$outlog
echo "" 
done


mkdir -p Dedup_idxstats

for file in $(ls Dedup/*.dedup.je.bam) 
do
echo "" 
echo "Processing file $file" 
base=$(basename "$file") 
output="Dedup_idxstats/$(echo ${base} |sed -e 's/.dedup.je.bam/.idxstats.txt/')" 
samtools index "$file" 
samtools idxstats $file > $output
done


multiqc -f \
Dedup_idxstats/*idxstats.txt \
-o ./ \
-n 07_Chromosome-mapping-post-dedup_MultiQC-report


outfile='Trimming/Demultiplex/Raw-read-counts.txt' 

for file in $(ls Trimming/Demultiplex/*.fastq.gz)
do
echo "File is $file"
base=$(basename "$file" .fastq.gz) 
count=$(echo $(zcat $file | wc -l) / 4 | bc) 
echo "Count is $count"
echo "$base $count">> $outfile
done

outfile='Trimming/Adapter/Trimmed-read-counts.txt' 

for file in $(ls Trimming/Adapter/*.fastq)
do
echo "File is $file"
base=$(basename "$file" _trimmed.fastq) 
count=$(echo $(expr $(cat $file | wc -l) / 4)) 
echo "Count is $count"
echo "$base $count">> $outfile
done


outfile='Alignment/Unique-map-read-counts.txt'

for dir in $(ls -d Alignment/*/) 
do
echo "" 
file=$dir/Aligned.sortedByCoord.out.bam 
base=$(basename "$dir") 
echo "File is $file" 
count=$(samtools view -c $file) 
echo "Count is $count"
echo "$base $count" >> $outfile
done


outfile='Dedup/Dedup-read-counts.txt'

for file in $(ls Dedup/*.dedup.je.bam) 
do
echo "" 
echo "File is $file" 
base=$(basename "$file" .dedup.je.bam) 
count=$(samtools view -c $file) 
echo "Count is $count"
echo "$base $count" >> $outfile
done


samtools faidx GENOMEDIR/HS.GRCh38.SINV.fa

cut -f1,2 GENOMEDIR/HS.GRCh38.SINV.fa.fai > GENOMEDIR/sizes.genome


mkdir -p Xlsite/collapsed
mkdir -p Xlsite/shifted/bed

for file in Dedup/*.dedup.je.bam 
do
echo "Processing $file" 
base=$(basename $file) 
outbed="Xlsite/collapsed/$(echo ${base} |sed -e 's/.dedup.je.bam/.collapsed.bed/')" 
outbedshifted="Xlsite/shifted/bed/$(echo ${base} |sed -e 's/.dedup.je.bam/.shifted.bed/')" 
bedtools bamtobed -i $file > $outbed 
bedtools shift -m 1 -p -1 -i $outbed -g GENOMEDIR/sizes.genome > $outbedshifted 
echo ""
done

mkdir -p Xlsite/shifted/bedgraph

for file in Xlsite/shifted/bed/*.shifted.bed 
do 
echo "Processing $file" 
base=$(basename $file .shifted.bed) 
outbgplus="Xlsite/shifted/bedgraph/$(echo ${base}).XL.RPM.+.bedgraph" 
outbgminus="Xlsite/shifted/bedgraph/$(echo ${base}).XL.RPM.-.bedgraph" 
TmpScale=$(bc <<< "scale=6;1000000/$(awk -v var="$base" '$1 == var { print $2}' Dedup/Dedup-read-counts.txt)") 
bedtools genomecov -bg -strand + -5 -i $file -g GENOMEDIR/sizes.genome -scale $TmpScale > $outbgplus 
bedtools genomecov -bg -strand - -5 -i $file -g GENOMEDIR/sizes.genome -scale $TmpScale > $outbgminus 
echo ""
done


mkdir -p Xlsite/shifted/bedgraph_sorted

LC_COLLATE=C

for file in Xlsite/shifted/bedgraph/*.bedgraph 
do 
echo "Processing $file" 
base=$(basename $file .bedgraph) 
output="Xlsite/shifted/bedgraph_sorted/$(echo ${base}).sorted.bedgraph" 
sort -k1,1 -k2,2n $file > $output
echo ""
done


mkdir -p Xlsite/shifted/bw

for file in Xlsite/shifted/bedgraph_sorted/*.bedgraph 
do 
echo "Processing $file" 
base=$(basename $file .sorted.bedgraph) 
output="Xlsite/shifted/bw/$(echo ${base}).bw" 
/software/kentUtils-v302.1.0/bin/linux.x86_64/bedGraphToBigWig –switches $file GENOMEDIR/sizes.genome $output
echo ""
done

mkdir -p DEW-seq/Extract

for file in Dedup/*.dedup.je.bam
do
echo "File is $file"
base=$(basename "$file")
output="DEW-seq/Extract/$(echo ${base} |sed -e 's/.dedup.je.bam/.extract.bed/')"
htseq-clip extract \
-i $file \
-e 1 \
-s s \
-g -1 \
--ignore \
-o $output
echo ""
done
 
mkdir -p DEW-seq/Count_w50s20

for file in DEW-seq/Extract/*.bed
do
echo "File is $file"
base=$(basename "$file")
output="DEW-seq/Count_w50s20/$(echo ${base} |sed -e 's/.extract.bed/.count.bed/')"
htseq-clip count \
-i $file \
-a GENOMEDIR/sliding_window/HS.GRCh38.SINV_attribute-fix.flattened.w50s20.txt.gz \
-o $output
echo ""
done


mkdir -p DEW-seq/Matrix_w50s20

htseq-clip createMatrix \
-i DEW-seq/Count_w50s20 \
-e .bed \
-o DEW-seq/Matrix_w50s20/Full_PURAB_iCLIP.matrix.txt.gz


