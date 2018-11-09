targetintrons="/home/bjf79/SpliceSeqPaperRawData/ScerSpliceSeq/BedFiles/TargetIntrons.sorted.bed"
transcriptsBed="/home/bjf79/SpliceSeqPaperRawData/ScerSpliceSeq/BedFiles/Supplemental_Table_S2_Booth_TSS.sorted.bed"
branchtothreess="/home/bjf79/SpliceSeqPaperRawData/ScerSpliceSeq/BedFiles/ScerBranchesTo3ss.sorted.bed" #the bed file needed to generate heatmap
SortedBamFile="/home/bjf79/SpliceSeqPaperRawData/ScerSpliceSeqFastq_ResultsOrig/STAR_DuplicateReadsRemovedAlignments/WTHansen/Filtered.bam"
ChromeSizesFile="/home/bjf79/SpliceSeqPaperRawData/ScerSpliceSeq/ScerReferenceGenome/SGD_Scer.fa.chrome.sizes"
R2Length=15

mkdir ./tempdir

#Count reads unambiguously lariat intermediates (R2 maps to branch window)
samtools view -b -f 131 -F 256 $SortedBamFile | bedtools coverage -c -b - -a <(bedtools shift -p $((5+R2Length)) -m $(($((5+R2Length)) * -1)) -i $branchtothreess -g $ChromeSizesFile | bedtools flank -s -l $((8 + R2Length)) -r 0 -i - -g $ChromeSizesFile) -F 1 -s -sorted -g $ChromeSizesFile > tempdir/UnambiguousLariatIntermediateCounts.bed

#Count reads that are ambiguous between pre-1st step and lariat-intermediate (R2 unmapped + R2 maps downstream of branch)
paste -d'\t' <(samtools view -b -f 131 -F 256 $SortedBamFile | bedtools coverage -c -b - -a <(bedtools slop -s -r $R2Length -l -5 -i $branchtothreess -g $ChromeSizesFile) -F 1 -s -sorted -g $ChromeSizesFile) <(samtools view -b -f 73 -F 256 $SortedBamFile | bedtools coverage -split -c -b - -a $branchtothreess -S -g $ChromeSizesFile -sorted) | awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7+$NF}' > tempdir/AmbiguousCounts.bed

#Count reads that are unambiguously pre-1st step (R1 is unspliced and R2 is neither ambiguous nor unambiguously lariat intermediate).
paste -d '\t' <(samtools view -b -f 65 -F 256 $SortedBamFile | bedtools coverage -split -c -b - -a $branchtothreess -S -g $ChromeSizesFile -sorted) tempdir/UnambiguousLariatIntermediateCounts.bed tempdir/AmbiguousCounts.bed | awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$14,$21,$7-$14-$21}' | bedtools intersect -loj -sorted -s -F 1 -b - -a <(samtools view -b -f 65 -F 256 $SortedBamFile | bedtools coverage -c -b - -a $targetintrons -S -g $ChromeSizesFile -sorted -split) | awk -F'\t' -v OFS='\t' '$8!="." {print $1,$2,$3,$4,$5,$6,$16,$14,$15} $8=="." {print $1,$2,$3,$4,$5,$6,"0","0",$7}' | sort -u -t$'\t' -k4,4 | bedtools sort -i - -faidx $ChromeSizesFile > EstimatedBranches.bed
