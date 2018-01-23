#!/usr/bin/env bash
#This script runs the STAR alignment before and after filtering reads for duplicates. It also makes a tab delimited file that counts spliced and unspliced read counts for each intron. Also plots some heatmaps

#==================
#Define relevant parameters as global variables for use throughout the script
#==================

genomefasta="$(dirname $0)/ScerReferenceGenome/SGD_Scer.fa" #if files are in home folder, do not use the tilde shortcut for your home folder. It will throw an error
genomegtf="$(dirname $0)/ScerReferenceGenome/SGD_Scer.gff" #Must be gff; if using gtf will have to change some parameters in the STAR genomeGenerate step
branchtothreess="$(dirname $0)/BedFiles/ScerBranchesTo3ss.sorted.bed" #the bed file needed to generate heatmap
targetintrons="$(dirname $0)/BedFiles/TargetIntrons.sorted.bed"
transcriptsBed="$(dirname $0)/BedFiles/Supplemental_Table_S2_Booth_TSS.sorted.bed"
fastqDirectory="$(dirname $0)/SampleData" #Directory should contain paired end read files, where each sample has an identical prefix and fastq files suffixed with _R1.fastq and _R2.fastq.

#==================
#Pre-alignment steps
#==================
set -xe #Debug mode

mkdir -p ./STAR_GenomeDirectory ./STAR_AllReadsAlignments ./STAR_DuplicateReadsRemovedAlignments ./Branchpoint_heatplots ./temp_files ./BedGraphCoverageFiles
samtools faidx $genomefasta
awk -v OFS='\t' '{print $1,$2}' $genomefasta.fai > $genomefasta.chrome.sizes
printf "Sample\tTargetReadCount\tOfftargetReadCount\tUnextendedPrimerReadCount\tUnmappedReadcount\tTargetUniqueReadCount\tOfftargetUniqueReadCount\tUnextendedPrimerUniqueReadCount\tUnmappedUniqueReadcount" > ./MappingStats.txt

#Make bed of all annotated intron-exon junctions (as 6 base intervals that flank the junction by three bases on each side, slop both -3 and flank both 6)
bedtools slop -b -3 -i $targetintrons -g $genomefasta.fai | bedtools flank -i stdin -b 6 -g $genomefasta.fai | sortBed -i stdin -faidx $genomefasta.fai > ./temp_files/intronexon_juncs.bed

#Make bed of transcription units that are of target genes
bedtools intersect -a $transcriptsBed -b $targetintrons -wa -sorted -g $genomefasta.chrome.sizes -s -F 1 > ./temp_files/TargetTranscripts.bed

#Make bed of regions where primers sit on target transcripts
bedtools flank -i $targetintrons -g $genomefasta.fai -l 0 -r 50 -s > ./PrimerRegions.sorted.bed

#Make genome files for STAR aligner
STAR --runMode genomeGenerate --genomeDir ./STAR_GenomeDirectory --genomeFastaFiles $genomefasta --sjdbGTFfile $genomegtf --sjdbOverhang 100 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS

#==================
#Alignment and processing loop
#==================

#Single pass mapping and additional analysis for each sample.
for myfilepath in $fastqDirectory/*_R1.fastq.gz
do
    samplename=$(basename -s "_R1.fastq.gz" $myfilepath)
    mkdir -p ./STAR_AllReadsAlignments/$samplename ./STAR_DuplicateReadsRemovedAlignments/$samplename
    echo $(date +"%b %d %T") .... intial STAR run of $samplename

    #Align reads
    STAR --genomeDir ./STAR_GenomeDirectory --readFilesCommand zcat --readFilesIn $fastqDirectory/${samplename}_R1.fastq $fastqDirectory/${samplename}_R2.fastq --readMapNumber -1 --alignIntronMin 20 --alignIntronMax 1100 --outFileNamePrefix ./STAR_AllReadsAlignments/$samplename/ --alignEndsType EndToEnd --clip3pAdapterSeq CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --clip5pNbases 7 0 --alignMatesGapMax 400 --alignSplicedMateMapLmin 16 --outSAMattributes All --runThreadN 4 --alignSJDBoverhangMin 1 --outSAMmultNmax 1 --outSAMunmapped Within KeepPairs --outFilterMismatchNmax 3
    samtools view -bS ./STAR_AllReadsAlignments/$samplename/Aligned.out.sam | samtools sort - > ./Aligned.out.bam
    samtools index ./Aligned.out.bam

    #Count target, nontarget, unextended primer, unmapped reads... Append to MappingStats.txt 
    MappedReadCount=$(samtools view -f 65 -F 260 Aligned.out.bam | wc -l)
    TotalPrimerRegionReadCount=$(samtools view -b -f 65 -F 256 Aligned.out.bam | bedtools coverage -b - -a PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S | awk -F'\t' '{sum+=$7} END{printf sum}')
    UnextendedPrimerReadCount=$(samtools view -b -f 65 -F 256 Aligned.out.bam | bedtools coverage -b - -a PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S -F 1 | awk -F'\t' '{sum+=$7} END{printf sum}')
    TargetReadCount=$(expr $TotalPrimerRegionReadCount - $UnextendedPrimerReadCount)
    OffTargetReadCount=$(expr $MappedReadCount - $TotalPrimerRegionReadCount)
    UnmappedReadCount=$(samtools view -f 69 -F 256 Aligned.out.bam | wc -l)
    printf "\n$samplename\t$TargetReadCount\t$OffTargetReadCount\t$UnextendedPrimerReadCount\t$UnmappedReadCount" >> ./MappingStats.txt
    
    #Filters out reads that aren't unique with respect to all bases in both reads in the original fastq files
    echo $(date +"%b %d %T") .... filtering PCR duplicate reads of $samplename
    paste <(zcat $fastqDirectory/${samplename}_R1.fastq | paste - - - -) <(zcat $fastqDirectory/${samplename}_R2.fastq | paste - - - -) | awk '{print $3$8, $0}' | sort -u -k1,1 - | awk -v OFS='\n' '{print $2,$4,$5,$6 > "./temp_files/MyR1s.filtered.fastq"; print $7,$9,$10,$11 > "./temp_files/MyR2s.filtered.fastq"}'

    #Align reads from the filtered bam file
    echo $(date +"%b %d %T") .... STAR run filtered reads of $samplename
    STAR --genomeDir ./STAR_GenomeDirectory --readFilesIn ./temp_files/MyR1s.filtered.fastq ./temp_files/MyR2s.filtered.fastq --alignIntronMin 20 --alignIntronMax 1100 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./STAR_DuplicateReadsRemovedAlignments/$samplename/ --alignEndsType EndToEnd --clip3pAdapterSeq CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --readMapNumber -1 --clip5pNbases 7 0 --alignMatesGapMax 400 --alignSplicedMateMapLmin 16 --outSAMattributes All --runThreadN 4 --alignSJDBoverhangMin 1 --outSAMmultNmax 1 --outSAMunmapped Within KeepPairs --outFilterMismatchNmax 3
    samtools index ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam

    #Count target, nontarget, unextended primer, unmapped reads... Append to MappingStats.txt
    MappedReadCount=$(samtools view -f 65 -F 260 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | wc -l)
    TotalPrimerRegionReadCount=$(samtools view -b -f 65 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools coverage -b - -a PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S | awk -F'\t' '{sum+=$7} END{printf sum}')
    UnextendedPrimerReadCount=$(samtools view -b -f 65 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools coverage -b - -a PrimerRegions.sorted.bed -sorted -g $genomefasta.chrome.sizes -S -F 1 | awk -F'\t' '{sum+=$7} END{printf sum}')
    TargetReadCount=$(expr $TotalPrimerRegionReadCount - $UnextendedPrimerReadCount)
    OffTargetReadCount=$(expr $MappedReadCount - $TotalPrimerRegionReadCount)
    UnmappedReadCount=$(samtools view -f 69 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | wc -l)
    printf "\t$TargetReadCount\t$OffTargetReadCount\t$UnextendedPrimerReadCount\t$UnmappedReadCount" >> MappingStats.txt
    
    #Make splicing-aware, paired-end gap-filled in bedgraph coverage files
    bedtools unionbedg -i <(samtools view -b -f 3 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -pc -strand -) <(samtools view -b -f 3 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -strand - -scale -1) <(samtools view -b -f 3 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -strand - -split) | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $4+$5+$6}' > ./STAR_DuplicateReadsRemovedAlignments/$samplename/pe_SpliceAwarePlusStrand.bedgraph
    bedtools unionbedg -i <(samtools view -b -f 3 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -pc -strand +) <(samtools view -b -f 3 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -strand + -scale -1) <(samtools view -b -f 3 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -strand + -split) | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, ($4+$5+$6)*-1}' > ./STAR_DuplicateReadsRemovedAlignments/$samplename/pe_SpliceAwareMinusStrand.bedgraph
    
    #Count IE (Intron-exon), EI, EE junctions for each intron and merge into single file
    samtools view -b -f 64 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools intersect -b - -a ./temp_files/intronexon_juncs.bed -S -sorted -g $genomefasta.chrome.sizes -split -wa -c -f 1 -bed | bedtools map -a $targetintrons -b stdin -s -c 7 -o mean -g $genomefasta.chrome.sizes -null 0 > ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronexonjuncCountsPerIntron.bed
    samtools view -b -f 64 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools coverage -b - -a $targetintrons -S -sorted -g $genomefasta.chrome.sizes -split > ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronCountsPerIntron.bed
    bedtools coverage -b ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam -a $transcriptsBed -s -sorted -g $genomefasta.chrome.sizes -split > ./STAR_DuplicateReadsRemovedAlignments/$samplename/CountsPerTranscript.bed
    paste -d "\t" ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronCountsPerIntron.bed ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronexonjuncCountsPerIntron.bed | awk -v OFS='\t' -F '\t' '{print $1,$2,$3,$4,$6,$7,$17}' | tee ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined.tab | awk -F '\t' -v OFS='\t' '{print $1"_"$2+1"_"$3"_"$5, $4, $6, $7}' | sort > ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined.temp.tab
    printf "chromosome\tstart\tstop\tstrand\tmotif\tannotated\tuniqe EE count\tmulti-mapping EE count\tEE-junction span\tJunction name\tIntronCount\t(EI+IE)/2\n" > ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined_SJout_merged.tab
    awk -v OFS='\t' -F '\t' '$4==2 {print $1"_"$2"_"$3"_-", $0} $4==1 {print $1"_"$2"_"$3"_+", $0}' ./STAR_DuplicateReadsRemovedAlignments/$samplename/SJ.out.tab | sort | join -a 1 -a 2 -e MISSING -t $'\t' -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.2,2.3,2.4 - ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined.temp.tab | sort -k1,1 -k2,2n >> ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined_SJout_merged.tab
    
    #Make bigwig coverage file for deepTools graphs
    samtools view -b -f 131 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools genomecov -ibam - -bg -5 > ./STAR_DuplicateReadsRemovedAlignments/$samplename/mybedgraph5.bedgraph
    /opt/kentUtils/bin/bedGraphToBigWig ./STAR_DuplicateReadsRemovedAlignments/$samplename/mybedgraph5.bedgraph $genomefasta.chrome.sizes ./STAR_DuplicateReadsRemovedAlignments/$samplename/ShortRead5ends.bw

    #Compute coverage for per-row normalization for deepTools graphs
    samtools view -b -f 67 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools coverage -b - -a $branchtothreess -sorted -g $genomefasta.chrome.sizes -split | awk -v OFS='\t' -F'\t' '$7>0 {print $1,$2,$3,$4,$7,$6}' > ./STAR_DuplicateReadsRemovedAlignments/$samplename/ScerBranchesTo3ss.counted.sorted.bed
    ##Filtered out transcripts for which the 5'ss is >400 bp away from TSS, since R2 reads for those genes wouldn't be expected to map
    samtools view -b -f 67 -F 256 ./STAR_DuplicateReadsRemovedAlignments/$samplename/Aligned.sortedByCoord.out.bam | bedtools coverage -b - -a <(bedtools flank -i $targetintrons -l 400 -r 0 -s -g $genomefasta.chrome.sizes | bedtools intersect -b - -a ./temp_files/TargetTranscripts.bed -g $genomefasta.chrome.sizes -wa -v -F 1) -sorted -S | awk -v OFS='\t' -F'\t' '$7>0 {print $1,$2,$3,$4,$7,$6}' > ./STAR_DuplicateReadsRemovedAlignments/$samplename/TargetTranscripts.counted.bed

    echo "...... making heat heatmap of 3'ends of reads around branches"

    #Branch Point relative CDF plot mean
    computeMatrix scale-regions -R ./STAR_DuplicateReadsRemovedAlignments/$samplename/ScerBranchesTo3ss.counted.sorted.bed -S ./STAR_DuplicateReadsRemovedAlignments/$samplename/ShortRead5ends.bw -m 20 --startLabel BranchA --endLabel 3SS -b 10 --unscaled5prime 10 -bs 1  -out ./STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz -p max/2 --missingDataAsZero --averageTypeBins sum -q
    introncount=$(zcat ./STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz | awk 'END {print NR}')
    zcat ./STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz | awk -F '\t' -v OFS='\t' 'NR==1 {print $0} NR>1 {sum=0; for(i=7; i<=NF; i++){sum+=$i}; s=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; runningsum=$5-sum; if(runningsum<0){runningsum=0}; for(i=7; i<=NF; i++){runningsum+=$i; s=s "\t" runningsum/$5*100}; print s}' | awk -F '\t' -v OFS='\t' 'NR==1 {print $0} NR>1 {print ($47-$7)*-1, $0}'| awk -F '\t' -v OFS='\t' 'NR==1{print $0;next}{print $0 | "sort -rn"}' | awk -F '\t' -v OFS='\t' 'NR==1{print $0;next}{$1=""; print $0}' | gzip - > ./STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrixNormalized.tab.gz
    plotHeatmap --heatmapHeight 4 --heatmapWidth 4 --yMax 100 --yMin 60 --colorMap Spectral --zMin 0 --zMax 100 --sortRegions descend --averageTypeSummaryPlot mean -m ./STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrixNormalized.tab.gz -out ./Branchpoint_heatplots/$samplename.BranchPointRelativeMeanCutoff.pdf --startLabel "A" --endLabel 3SS --x "cDNA Relative End Locations" -y "Relative % cDNA ends" --zMax 100 --regionsLabel "$introncount transcripts"
    gunzip -f ./STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrixNormalized.tab.gz

    #Branch Point absolute PDF plot mean
    plotHeatmap --heatmapHeight 4 --heatmapWidth 4 --colorMap Reds --yMin 0 --yMax 30 --sortRegions descend --averageTypeSummaryPlot mean -m ./STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz -out ./Branchpoint_heatplots/$samplename.BranchPointAbsoluteMeanPeak.pdf --startLabel "A" --endLabel 3SS --x "cDNA End Locations" -y "Meta-region mean" --regionsLabel "$introncount transcripts"
    gunzip -f ./STAR_DuplicateReadsRemovedAlignments/$samplename/BranchPointCountMatrix.tab.gz

    #TSS relative CDF plot mean
    computeMatrix reference-point -R ./STAR_DuplicateReadsRemovedAlignments/$samplename/TargetTranscripts.counted.bed -S ./STAR_DuplicateReadsRemovedAlignments/$samplename/ShortRead5ends.bw -a 500 -b 500 -bs 1 -out ./STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz -p max/2 --missingDataAsZero -q
    introncount=$(zcat ./STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz | awk 'END {print NR}')
    zcat ./STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz | awk -F '\t' -v OFS='\t' 'NR==1 {print $0} NR>1 {sum=0; for(i=7; i<=NF; i++){sum+=$i}; s=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; runningsum=$5-sum; if(runningsum<0){runningsum=0}; for(i=7; i<=NF; i++){runningsum+=$i; s=s "\t" runningsum/$5*100}; print s}' | awk -F '\t' -v OFS='\t' 'NR==1 {print $0} NR>1 {print ($47-$7)*-1, $0}'| awk -F '\t' -v OFS='\t' 'NR==1{print $0;next}{print $0 | "sort -rn"}' | awk -F '\t' -v OFS='\t' 'NR==1{print $0;next}{$1=""; print $0}' | gzip - > ./STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrixNormalized.tab.gz
    plotHeatmap --heatmapHeight 4 --heatmapWidth 4 --yMax 100 --yMin 0 --colorMap jet --zMin 0 --zMax 100 --sortRegions descend --averageTypeSummaryPlot mean -m ./STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrixNormalized.tab.gz -out ./Branchpoint_heatplots/$samplename.TSSRelativeMeanCutoff.pdf --startLabel "TSS" --x "cDNA Relative End Locations" -y "Relative % cDNA ends" --zMax 100 --regionsLabel "$introncount transcripts"
    gunzip -f ./STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrixNormalized.tab.gz

    #TSS absolute PDF plot mean
    plotHeatmap --heatmapHeight 4 --heatmapWidth 4 --zMin 0 --zMax 100 --yMax 1500 --yMin 0 --sortUsing sum --colorMap Blues --sortRegions descend --averageTypeSummaryPlot mean -m ./STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz -out ./Branchpoint_heatplots/$samplename.TSSAbsoluteMeanPeak.pdf --startLabel "TSS" --x "cDNA End Locations" -y "Meta-region mean" --regionsLabel "$introncount transcripts"
    gunzip -f ./STAR_DuplicateReadsRemovedAlignments/$samplename/TSSCountMatrix.tab.gz 

rm ./Aligned.out.bam ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronCountsPerIntron.bed ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronexonjuncCountsPerIntron.bed ./STAR_AllReadsAlignments/$samplename/Aligned.out.sam ./STAR_DuplicateReadsRemovedAlignments/$samplename/intronCoverageCombined.temp.tab
done

rm -R ./temp_files
