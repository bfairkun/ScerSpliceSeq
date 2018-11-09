#!/usr/bin/env bash
#This script runs the STAR alignment before and after filtering reads for duplicates. It also makes a tab delimited file that counts spliced and unspliced read counts for each intron. Also plots some heatmaps

#==================
#Define relevant parameters as global variables for use throughout the script
#==================

R2Length=15 #Length of Read2. relevant to determine window size for counting reads that lie within a region indicative of a branched lariat intermediate
genomefasta="$(dirname $0)/ScerReferenceGenome/SGD_Scer.fa" #if files are in home folder, do not use the tilde shortcut for your home folder. It will throw an error
genomegtf="$(dirname $0)/ScerReferenceGenome/SGD_Scer.gff" #Must be gff; if using gtf will have to change some parameters in the STAR genomeGenerate step
branchtothreess="$(dirname $0)/BedFiles/ScerBranchesTo3ss.filtered.bed" #the bed file needed to generate heatmap
targetintrons="$(dirname $0)/BedFiles/TargetIntrons.sorted.bed"
transcriptsBed="$(dirname $0)/BedFiles/Supplemental_Table_S2_Booth_TSS.sorted.bed"
fastqDirectory="./MikesFastq" #Directory should contain paired end read files, where each sample has an identical prefix and fastq files suffixed with _R1.fastq and _R2.fastq.
OutDir="./MikesFastq_Results"

#==================
#Pre-alignment steps
#==================
set -xe #Debug mode

mkdir -p $OutDir $OutDir/STAR_GenomeDirectory $OutDir/STAR_AllReadsAlignments $OutDir/STAR_DuplicateReadsRemovedAlignments $OutDir/Branchpoint_heatplots $OutDir/temp_files $OutDir/BedGraphCoverageFiles
samtools faidx $genomefasta
awk -v OFS='\t' '{print $1,$2}' $genomefasta.fai > $genomefasta.chrome.sizes
printf "Sample\tTargetReadCount\tOfftargetReadCount\tUnextendedPrimerReadCount\tUnmappedReadcount\tTargetUniqueReadCount\tOfftargetUniqueReadCount\tUnextendedPrimerUniqueReadCount\tUnmappedUniqueReadcount" > $OutDir/MappingStats.txt

#Make bed of all annotated intron-exon junctions (as 6 base intervals that flank the junction by three bases on each side, slop both -3 and flank both 6)
bedtools slop -b -3 -i $targetintrons -g $genomefasta.fai | bedtools flank -i stdin -b 6 -g $genomefasta.fai | sortBed -i stdin -faidx $genomefasta.fai > $OutDir/temp_files/intronexon_juncs.bed

#Make bed of transcription units that are of target genes
bedtools intersect -a $transcriptsBed -b $targetintrons -wa -sorted -g $genomefasta.chrome.sizes -s -F 1 > $OutDir/temp_files/TargetTranscripts.bed

#Make bed of regions where primers sit on target transcripts
bedtools flank -i $targetintrons -g $genomefasta.fai -l 0 -r 50 -s > $OutDir/PrimerRegions.sorted.bed

#Make genome files for STAR aligner
STAR --runMode genomeGenerate --genomeDir $OutDir/STAR_GenomeDirectory --genomeFastaFiles $genomefasta --sjdbGTFfile $genomegtf --sjdbOverhang 100 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS
#==================
#Alignment and processing loop
#==================

#Single pass mapping and additional analysis for each sample.
for myfilepath in $fastqDirectory/*_R1.fastq.gz
do
    samplename=$(basename -s "_R1.fastq.gz" $myfilepath)
    mkdir -p $OutDir/STAR_AllReadsAlignments/$samplename $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename
    echo $(date +"%b %d %T") .... intial STAR run of $samplename

    #Align reads
    STAR --genomeDir $OutDir/STAR_GenomeDirectory --readFilesCommand zcat --readFilesIn $fastqDirectory/${samplename}_R1.fastq.gz $fastqDirectory/${samplename}_R2.fastq.gz --readMapNumber -1 --alignIntronMin 20 --alignIntronMax 1100 --outFileNamePrefix $OutDir/STAR_AllReadsAlignments/$samplename/ --alignEndsType EndToEnd --clip3pAdapterSeq CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --clip5pNbases 7 0 --alignMatesGapMax 400 --alignSplicedMateMapLmin 16 --outSAMattributes All --runThreadN 4 --alignSJDBoverhangMin 1 --outSAMmultNmax 1 --outSAMunmapped Within KeepPairs --outFilterMismatchNmax 3
    samtools view -bS $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.out.sam | samtools sort - > $OutDir/Aligned.out.bam
    samtools index $OutDir/Aligned.out.bam

    #Count target, nontarget, unextended primer, unmapped reads... Append to MappingStats.txt
    MappedReadCount=$(samtools view -f 65 -F 260 $OutDir/Aligned.out.bam | wc -l)
    TotalPrimerRegionReadCount=$(samtools view -b -f 65 -F 260 $OutDir/Aligned.out.bam | bedtools intersect -abam - -wa -b $OutDir/PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S -ubam | samtools view | wc -l)
    UnextendedPrimerReadCount=$(samtools view -b -f 65 -F 260 $OutDir/Aligned.out.bam | bedtools intersect -abam - -wa -b $OutDir/PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S -ubam -f 1 | samtools view | wc -l)
    TargetReadCount=$(expr $TotalPrimerRegionReadCount - $UnextendedPrimerReadCount)
    OffTargetReadCount=$(expr $MappedReadCount - $TotalPrimerRegionReadCount)
    OfftargetReadCountFromUnmapped=$(STAR --genomeDir $OutDir/STAR_GenomeDirectory --readFilesIn <(samtools bam2fq -f 68 $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.out.sam) <(samtools bam2fq -f 136 $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.out.sam) --readMapNumber -1 --alignIntronMin 20 --alignIntronMax 1100 --outFilterMultimapNmax 1000 --alignEndsType EndToEnd --clip3pAdapterSeq CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --clip5pNbases 32 0 --alignMatesGapMax 400 --alignSplicedMateMapLmin 16 --outSAMattributes All --runThreadN 4 --alignSJDBoverhangMin 1 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outSAMunmapped Within | tee SecondPass.bam | samtools view -b -f 65 -F 260 | bedtools intersect -abam - -wa -b $OutDir/PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S -ubam -v | samtools view | wc -l)
    TotalOffTargetReadCount=$(expr $OffTargetReadCount + $OfftargetReadCountFromUnmapped)
    UnmappedReadCount=$(expr $(samtools view -f 69 -F 256 $OutDir/Aligned.out.bam | wc -l) - $OfftargetReadCountFromUnmapped)
    printf "\n$samplename\t$TargetReadCount\t$TotalOffTargetReadCount\t$UnextendedPrimerReadCount\t$UnmappedReadCount" >> $OutDir/MappingStats.txt

    #Filters out reads that aren't unique with respect to all bases in both reads in the original fastq files
    echo $(date +"%b %d %T") .... filtering PCR duplicate reads of $samplename
    paste <(zcat $fastqDirectory/${samplename}_R1.fastq | paste - - - -) <(zcat $fastqDirectory/${samplename}_R2.fastq | paste - - - -) | awk '{print $3$8, $0}' | sort -u -k1,1 - | awk -v OutDir="$OutDir" -v OFS='\n' '{print $2,$4,$5,$6 > OutDir"/temp_files/MyR1s.filtered.fastq"; print $7,$9,$10,$11 > OutDir"/temp_files/MyR2s.filtered.fastq"}'

    #Align reads from the filtered bam file
    echo $(date +"%b %d %T") .... STAR run filtered reads of $samplename
    STAR --genomeDir $OutDir/STAR_GenomeDirectory --readFilesIn $OutDir/temp_files/MyR1s.filtered.fastq $OutDir/temp_files/MyR2s.filtered.fastq --alignIntronMin 20 --alignIntronMax 1100 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/ --alignEndsType EndToEnd --clip3pAdapterSeq CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --readMapNumber -1 --clip5pNbases 7 0 --alignMatesGapMax 400 --alignSplicedMateMapLmin 16 --outSAMattributes All --runThreadN 4 --alignSJDBoverhangMin 1 --outSAMmultNmax 1 --outSAMunmapped Within KeepPairs --outFilterMismatchNmax 3
    samtools index $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam

    #Count target, nontarget, unextended primer, unmapped reads... Append to MappingStats.txt
    MappedReadCount=$(samtools view -f 65 -F 260 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | wc -l)
    TotalPrimerRegionReadCount=$(samtools view -b -f 65 -F 260 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools intersect -abam - -wa -b $OutDir/PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S -ubam | samtools view | wc -l)
    UnextendedPrimerReadCount=$(samtools view -b -f 65 -F 260 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools intersect -abam - -wa -b $OutDir/PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S -ubam -f 1 | samtools view | wc -l)
    TargetReadCount=$(expr $TotalPrimerRegionReadCount - $UnextendedPrimerReadCount)
    OffTargetReadCount=$(expr $MappedReadCount - $TotalPrimerRegionReadCount)
    OfftargetReadCountFromUnmapped=$(STAR --genomeDir $OutDir/STAR_GenomeDirectory --readFilesIn <(samtools bam2fq -f 68 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam) <(samtools bam2fq -f 136 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam) --readMapNumber -1 --alignIntronMin 20 --alignIntronMax 1100 --alignEndsType EndToEnd --clip3pAdapterSeq CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --clip5pNbases 32 0 --alignMatesGapMax 400 --alignSplicedMateMapLmin 16 --outSAMattributes All --runThreadN 4 --alignSJDBoverhangMin 1 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate | samtools view -b -f 65 -F 260 | bedtools intersect -abam - -wa -b $OutDir/PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S -ubam -v | samtools view | wc -l)
    TotalOffTargetReadCount=$(expr $OffTargetReadCount + $OfftargetReadCountFromUnmapped)
    UnmappedReadCount=$(expr $(samtools view -f 69 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | wc -l) - $OfftargetReadCountFromUnmapped)
    printf "\t$TargetReadCount\t$TotalOffTargetReadCount\t$UnextendedPrimerReadCount\t$UnmappedReadCount" >> $OutDir/MappingStats.txt
    
    #Make splicing-aware, paired-end gap-filled in bedgraph coverage files
    bedtools unionbedg -i <(samtools view -b -f 3 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -pc -strand -) <(samtools view -b -f 3 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -strand - -scale -1) <(samtools view -b -f 3 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -strand - -split) | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $4+$5+$6}' > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/pe_SpliceAwarePlusStrand.bedgraph
    bedtools unionbedg -i <(samtools view -b -f 3 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -pc -strand +) <(samtools view -b -f 3 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -strand + -scale -1) <(samtools view -b -f 3 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -strand + -split) | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, ($4+$5+$6)*-1}' > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/pe_SpliceAwareMinusStrand.bedgraph

    #Filter out alignments that are under 30bp aligned (i.e. either observed template length (SAM field 9), or for unpaired alignments the CIGAR matches are > 30
    cat <(samtools view -H $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam) <(samtools view $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | perl -F'\t' -lane 'BEGIN {use List::Util qw(sum);} @mymatches = $F[5] =~ m/(\d+)M/g; if( (abs($F[8]) >= 30) || (sum(@mymatches)>30) ) {print}') | samtools view -bS - > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam
    samtools index $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam

    #Count IE (Intron-exon), EI, EE junctions for each intron and merge into single file. Also count pre1st step unspliced, post1st step unspliced, and ambiguous unspliced for introns with well defined branches
    samtools view -b -f 64 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam | bedtools intersect -b - -a $OutDir/temp_files/intronexon_juncs.bed -S -sorted -g $genomefasta.chrome.sizes -split -wa -c -f 1 -bed | bedtools map -a $targetintrons -b stdin -s -c 7 -o mean -g $genomefasta.chrome.sizes -null 0 > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronexonjuncCountsPerIntron.bed
    samtools view -b -f 131 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam | bedtools coverage -c -b - -a <(bedtools shift -p $((5+R2Length)) -m $(($((5+R2Length)) * -1)) -i $branchtothreess -g $genomefasta.chrome.sizes | bedtools flank -s -l $((8 + R2Length)) -r 0 -i - -g $genomefasta.chrome.sizes) -F 1 -s -sorted -g $genomefasta.chrome.sizes > $OutDir/temp_files/UnambiguousLariatIntermediateCounts.bed
    paste -d'\t' <(samtools view -b -f 131 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam | bedtools coverage -c -b - -a <(bedtools slop -s -r $R2Length -l -5 -i $branchtothreess -g $genomefasta.chrome.sizes) -F 1 -s -sorted -g $genomefasta.chrome.sizes) <(samtools view -b -f 73 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam | bedtools coverage -split -c -b - -a $branchtothreess -S -g $genomefasta.chrome.sizes -sorted) | awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7+$NF}' > $OutDir/temp_files/AmbiguousCounts.bed
    paste -d '\t' <(samtools view -b -f 65 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam | bedtools coverage -split -c -b - -a $branchtothreess -S -g $genomefasta.chrome.sizes -sorted) $OutDir/temp_files/UnambiguousLariatIntermediateCounts.bed $OutDir/temp_files/AmbiguousCounts.bed | tee mytee.txt | awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$14,$21,$7-$14-$21}' | bedtools intersect -loj -sorted -s -F 1 -b - -a <(samtools view -b -f 65 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam | bedtools coverage -split -c -b - -a $targetintrons -S -g $genomefasta.chrome.sizes -sorted) | awk -F'\t' -v OFS='\t' '$8!="." {print $1,$2,$3,$4,$5,$6,$16,$14,$15} $8=="." {print $1,$2,$3,$4,$5,$6,"0","0",$7}' | sort -u -t$'\t' -k4,4 | bedtools sort -i - -faidx $genomefasta.chrome.sizes > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronCountsPerIntron.bed
    bedtools coverage -b $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam -a $transcriptsBed -s -sorted -g $genomefasta.chrome.sizes -split > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/CountsPerTranscript.bed
    paste -d '\t' $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronCountsPerIntron.bed $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronexonjuncCountsPerIntron.bed | awk -v OFS='\t' -F '\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$19}' | tee $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined.tab | awk -F '\t' -v OFS='\t' '{print $1"_"$2+1"_"$3"_"$6, $4, $7, $8, $9, $10}' | sort > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined.temp.tab
    printf "chromosome\tstart\tstop\tstrand\tmotif\tannotated\tuniqe EE count\tmulti-mapping EE count\tEE-junction span\tJunction name\tPreFirstStep\tLariatIntermediate\tAmbiguousUnspliced\t(EI+IE)/2\n" > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined_SJout_merged.tab
    awk -v OFS='\t' -F '\t' '$4==2 {print $1"_"$2"_"$3"_-", $0} $4==1 {print $1"_"$2"_"$3"_+", $0}' $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/SJ.out.tab | sort | join -a 1 -a 2 -e MISSING -t $'\t' -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.2,2.3,2.4,2.5,2.6 - $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined.temp.tab | sort -k1,1 -k2,2n >> $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined_SJout_merged.tab

    #Make bigwig coverage file for deepTools graphs
    samtools view -b -f 131 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam | bedtools genomecov -ibam - -bg -5 > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/mybedgraph5.bedgraph
    /opt/kentUtils/bin/bedGraphToBigWig $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/mybedgraph5.bedgraph $genomefasta.chrome.sizes $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/ShortRead5ends.bw

    #Compute coverage for per-row normalization for deepTools graphs
    samtools view -b -f 67 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam | bedtools coverage -b - -a $branchtothreess -sorted -g $genomefasta.chrome.sizes -split | awk -v OFS='\t' -F'\t' '$7>0 {print $1,$2,$3,$4,$7,$6}' > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/ScerBranchesTo3ss.counted.sorted.bed
    ##Filtered out transcripts for which the 5'ss is >400 bp away from TSS, since R2 reads for those genes wouldn't be expected to map
    samtools view -b -f 67 -F 256 $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/Filtered.bam | bedtools coverage -b - -a <(bedtools flank -i $targetintrons -l 400 -r 0 -s -g $genomefasta.chrome.sizes | bedtools intersect -b - -a $OutDir/temp_files/TargetTranscripts.bed -g $genomefasta.chrome.sizes -wa -v -F 1) -sorted -S | awk -v OFS='\t' -F'\t' '$7>0 {print $1,$2,$3,$4,$7,$6}' > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TargetTranscripts.counted.bed

    echo "...... making heat heatmap of 3'ends of reads around branches"

    #Branch Point relative CDF plot mean
    computeMatrix scale-regions -R $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/ScerBranchesTo3ss.counted.sorted.bed -S $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/ShortRead5ends.bw -m 20 --startLabel BranchA --endLabel 3SS -b 10 --unscaled5prime 10 -bs 1  -out $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz -p max/2 --missingDataAsZero --averageTypeBins sum -q
    introncount=$(zcat $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz | awk 'END {print NR}')
    zcat $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz | awk -F '\t' -v OFS='\t' 'NR==1 {print $0} NR>1 {sum=0; for(i=7; i<=NF; i++){sum+=$i}; s=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; runningsum=$5-sum; if(runningsum<0){runningsum=0}; for(i=7; i<=NF; i++){runningsum+=$i; s=s "\t" runningsum/$5*100}; print s}' | awk -F '\t' -v OFS='\t' 'NR==1 {print $0} NR>1 {print ($47-$7)*-1, $0}'| awk -F '\t' -v OFS='\t' 'NR==1{print $0;next}{print $0 | "sort -rn"}' | awk -F '\t' -v OFS='\t' 'NR==1{print $0;next}{$1=""; print $0}' | gzip - > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrixNormalized.tab.gz
    plotHeatmap --heatmapHeight 4 --heatmapWidth 4 --yMax 100 --yMin 60 --colorMap Spectral --zMin 0 --zMax 100 --sortRegions descend --averageTypeSummaryPlot median -m $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrixNormalized.tab.gz -out $OutDir/Branchpoint_heatplots/$samplename.BranchPointRelativeMeanCutoff.pdf --startLabel "A" --endLabel 3SS --x "intronic coverage (relative)" -y "Meta-region mean" --zMax 100 --regionsLabel "$introncount transcripts"
    gunzip -f $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrixNormalized.tab.gz

    #Branch Point absolute PDF plot mean
    plotHeatmap --heatmapHeight 4 --heatmapWidth 4 --colorMap Reds --yMin 0 --yMax 30 --sortRegions descend --averageTypeSummaryPlot median -m $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz -out $OutDir/Branchpoint_heatplots/$samplename.BranchPointAbsoluteMeanPeak.pdf --startLabel "A" --endLabel 3SS --x "coverage" -y "Meta-region mean" --regionsLabel "$introncount transcripts"
    gunzip -f $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz

    #TSS relative CDF plot mean
    computeMatrix reference-point -R $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TargetTranscripts.counted.bed -S $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/ShortRead5ends.bw -a 500 -b 500 -bs 1 -out $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz -p max/2 --missingDataAsZero -q
    introncount=$(zcat $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz | awk 'END {print NR}')
    zcat $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz | awk -F '\t' -v OFS='\t' 'NR==1 {print $0} NR>1 {sum=0; for(i=7; i<=NF; i++){sum+=$i}; s=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; runningsum=$5-sum; if(runningsum<0){runningsum=0}; for(i=7; i<=NF; i++){runningsum+=$i; s=s "\t" runningsum/$5*100}; print s}' | awk -F '\t' -v OFS='\t' 'NR==1 {print $0} NR>1 {print ($47-$7)*-1, $0}'| awk -F '\t' -v OFS='\t' 'NR==1{print $0;next}{print $0 | "sort -rn"}' | awk -F '\t' -v OFS='\t' 'NR==1{print $0;next}{$1=""; print $0}' | gzip - > $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrixNormalized.tab.gz
    plotHeatmap --heatmapHeight 4 --heatmapWidth 4 --yMax 100 --yMin 0 --colorMap jet --zMin 0 --zMax 100 --sortRegions descend --averageTypeSummaryPlot mean -m $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrixNormalized.tab.gz -out $OutDir/Branchpoint_heatplots/$samplename.TSSRelativeMeanCutoff.pdf --startLabel "TSS" --x "transcript coverage (relative)" -y "Meta-region mean" --zMax 100 --regionsLabel "$introncount transcripts"
    gunzip -f $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrixNormalized.tab.gz

    #TSS absolute PDF plot mean
    plotHeatmap --heatmapHeight 4 --heatmapWidth 4 --zMin 0 --zMax 100 --yMax 1500 --yMin 0 --sortUsing sum --colorMap Blues --sortRegions descend --averageTypeSummaryPlot mean -m $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz -out $OutDir/Branchpoint_heatplots/$samplename.TSSAbsoluteMeanPeak.pdf --startLabel "TSS" --x "coverage" -y "Meta-region mean" --regionsLabel "$introncount transcripts"
    gunzip -f $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz
    rm $OutDir/Aligned.out.bam
    rm $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronCountsPerIntron.bed $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronexonjuncCountsPerIntron.bed $OutDir/STAR_AllReadsAlignments/$samplename/Aligned.out.sam $OutDir/STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined.temp.tab
done

rm -R $OutDir/temp_files
