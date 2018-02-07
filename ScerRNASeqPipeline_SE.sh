#!/usr/bin/env bash
#This script runs the STAR alignment before and after filtering reads for duplicates. It also makes a tab delimited file that counts spliced and unspliced read counts for each intron. Also plots some heatmaps

#==================
#Define relevant parameters as global variables for use throughout the script
#==================

genomefasta="$(dirname $0)/ScerReferenceGenome/SGD_Scer.fa" #if files are in home folder, do not use the tilde shortcut for your home folder. It will throw an error
genomegtf="$(dirname $0)/ScerReferenceGenome/SGD_Scer.gff" #Must be gff; if using gtf will have to change some parameters in the STAR genomeGenerate step
targetintrons="$(dirname $0)/BedFiles/TargetIntrons.sorted.bed"
transcriptsBed="$(dirname $0)/BedFiles/Supplemental_Table_S2_Booth_TSS.sorted.bed"
fastqDirectory="./ScerRNAseqFastq" #Directory should contain paired end read files, where each sample has an identical prefix and fastq files suffixed with _R1.fastq and _R2.fastq.
OutDir="./ScerRNAseqFastq_Results"

#==================
#Pre-alignment steps
#==================
set -xe #Debug mode

mkdir -p $OutDir $OutDir/STAR_GenomeDirectory $OutDir/STAR_AllReadsAlignments $OutDir/temp_files
awk -v OFS='\t' '{print $1,$2}' $genomefasta.fai > $genomefasta.chrome.sizes
printf "Sample\tTargetReadCount\tOfftargetReadCount\tUnextendedPrimerReadCount\tUnmappedReadcount" > $OutDir/MappingStats.txt

#Make bed of all annotated intron-exon junctions (as 6 base intervals that flank the junction by three bases on each side, slop both -3 and flank both 6)
bedtools slop -b -3 -i $targetintrons -g $genomefasta.fai | bedtools flank -i stdin -b 6 -g $genomefasta.fai | sortBed -i stdin -faidx $genomefasta.fai > $OutDir/temp_files/intronexon_juncs.bed

#Make bed of transcription units that are of target genes
bedtools intersect -a $transcriptsBed -b $targetintrons -wa -sorted -g $genomefasta.chrome.sizes -s -F 1 > $OutDir/temp_files/TargetTranscripts.bed

#Make genome files for STAR aligner
STAR --runMode genomeGenerate --genomeDir $OutDir/STAR_GenomeDirectory --genomeFastaFiles $genomefasta --sjdbGTFfile $genomegtf --sjdbOverhang 100 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS

#==================
#Alignment and processing loop
#==================

#Single pass mapping and additional analysis for each sample.
for myfilepath in $fastqDirectory/*fastq.gz
do
    samplename=$(basename -s ".fastq.gz" $myfilepath)
    mkdir -p $OutDir/STAR_AllReadsAlignments/$samplename
    echo $(date +"%b %d %T") .... intial STAR run of $samplename
    ReadLength=$(zcat $fastqDirectory/${samplename}.fastq.gz | head -2 | awk 'NR==2 {print length($1)}')

    #Align reads
    STAR --genomeDir $OutDir/STAR_GenomeDirectory --readFilesCommand zcat --readFilesIn $fastqDirectory/${samplename}.fastq.gz --readMapNumber -1 --alignIntronMin 20 --alignIntronMax 1100 --outFileNamePrefix $OutDir/STAR_AllReadsAlignments/$samplename/ --alignEndsType EndToEnd --clip3pAdapterSeq CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --clip5pNbases 0 --alignMatesGapMax 400 --alignSplicedMateMapLmin 16 --outSAMattributes All --runThreadN 4 --alignSJDBoverhangMin 1 --outSAMmultNmax 1 --outSAMunmapped Within --outFilterMismatchNmax 3 --outWigType bedGraph --outSAMtype BAM SortedByCoordinate --outWigNorm None
    samtools index $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.sortedByCoord.out.bam

    #Count target, nontarget, unextended primer, unmapped reads... Append to MappingStats.txt
    MappedReadCount=$(samtools view -F 260 $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.sortedByCoord.out.bam | wc -l)
    TargetReadCount=$(samtools view -b -F 256 $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools intersect -abam - -b <(bedtools slop -i $targetintrons -b 3 -g $genomefasta.fai) -sorted -g $genomefasta.chrome.sizes -S -ubam | samtools view | wc -l)
    OffTargetReadCount=$(expr $MappedReadCount - $TargetReadCount)
    UnmappedReadCount=$(samtools view -f 4 -F 256 $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.sortedByCoord.out.bam | wc -l)
    printf "\n$samplename\t$TargetReadCount\t$OffTargetReadCount\tNA\t$UnmappedReadCount" >> $OutDir/MappingStats.txt

    #Count IE (Intron-exon), EI, EE junctions for each intron and merge into single file
    samtools view -b -F 256 $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools intersect -b - -a $OutDir/temp_files/intronexon_juncs.bed -S -sorted -g $genomefasta.chrome.sizes -split -wa -c -f 1 -bed | bedtools map -a $targetintrons -b stdin -s -c 7 -o mean -g $genomefasta.chrome.sizes -null 0 > $OutDir/STAR_AllReadsAlignments/$samplename/intronexonjuncCountsPerIntron.bed
    samtools view -b -F 256 $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools coverage -b - -a <(bedtools slop -i $targetintrons -b -3 -g $genomefasta.fai) -S -sorted -g $genomefasta.chrome.sizes -split | bedtools slop -i - -b 3 -g $genomefasta.fai > $OutDir/STAR_AllReadsAlignments/$samplename/intronCountsPerIntron.bed
    bedtools coverage -b $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.sortedByCoord.out.bam -a $transcriptsBed -s -sorted -g $genomefasta.chrome.sizes -split > $OutDir/STAR_AllReadsAlignments/$samplename/CountsPerTranscript.bed
    paste -d "\t" $OutDir/STAR_AllReadsAlignments/$samplename/intronCountsPerIntron.bed $OutDir/STAR_AllReadsAlignments/$samplename/intronexonjuncCountsPerIntron.bed | awk -v OFS='\t' -F '\t' '{print $1,$2,$3,$4,$6,$7,$17}' | tee $OutDir/STAR_AllReadsAlignments/$samplename/intronCoverageCombined.tab | awk -F '\t' -v OFS='\t' '{print $1"_"$2+1"_"$3"_"$5, $4, $6, $7}' | sort > $OutDir/STAR_AllReadsAlignments/$samplename/intronCoverageCombined.temp.tab
    printf "chromosome\tstart\tstop\tstrand\tmotif\tannotated\tuniqe EE count\tmulti-mapping EE count\tEE-junction span\tJunction name\tIntronCount\t(EI+IE)/2\tIntronCountsPerKbMappingSpace\tExonExonJunctionCountsPerKbMappingSpace\tPercentUnspliced\n" > $OutDir/STAR_AllReadsAlignments/$samplename/intronCoverageCombined_SJout_merged.tab
    awk -v OFS='\t' -F '\t' '$4==2 {print $1"_"$2"_"$3"_-", $0} $4==1 {print $1"_"$2"_"$3"_+", $0}' $OutDir/STAR_AllReadsAlignments/$samplename/SJ.out.tab | sort | join -a 1 -a 2 -e MISSING -t $'\t' -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.2,2.3,2.4 - $OutDir/STAR_AllReadsAlignments/$samplename/intronCoverageCombined.temp.tab | sort -k1,1 -k2,2n | awk -F'\t' -v ReadLength="$ReadLength" -v OFS='\t' '{Spliced=$7/(2*(ReadLength-3))*1000; Unspliced=$11/(2*(ReadLength-3)+$3-$2)*1000; print $0, Unspliced, Spliced, Unspliced/(Spliced+Unspliced)}' >> $OutDir/STAR_AllReadsAlignments/$samplename/intronCoverageCombined_SJout_merged.tab

    rm $OutDir/STAR_AllReadsAlignments/$samplename/intronCountsPerIntron.bed $OutDir/STAR_AllReadsAlignments/$samplename/intronexonjuncCountsPerIntron.bed $OutDir/STAR_AllReadsAlignments/$samplename/intronCoverageCombined.temp.tab
done

rm -R $OutDir/temp_files
