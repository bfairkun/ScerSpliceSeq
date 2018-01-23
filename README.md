# Screen-Seq Analysis sample tutorial

The bash script in this repository serves as a pipeline to analyze paired-end (PE) splice-seq data S. cerevisiae. Modify the script as needed. Will add more description to this readme later.

## Dependencies

The following softwares should be available from PATH environment
- [STAR](https://github.com/alexdobin/STAR)
- [samtools](http://samtools.sourceforge.net)
- [bedtools](http://bedtools.readthedocs.io/en/latest/#)
- [deepTools](http://deeptools.readthedocs.io/en/latest/)

## Editing and running the script

Throughout the bash script (script_master_copy-3.sh), output directories and files will be generated. In its current form, all output directories will be made from from the directory from which the script is called. Therefore, I recommend making a new directory to house all of the output, and then cloning this repository within that directory, and running it **from** that directory that you made...

```bash
mkdir MyNewDirectory
cd MyNewDirectory
git clone myrepo
```

But before running the script, you should probably check the global variables at the top of the script that define some filepaths. In its current form, all of these are defined as paths to files that come with this repository, so the script should run fine as is. However, if you want to run the script on your own fastq fata (as opposed to the test data that comes with this repository), you will need to change the fastq directory filepath. Note that the fastq files in your fastq directory must all be paired, compressed, and named as {samplename}\_R1.fastq.gz and {samplename}_R2.fastq.gz. Otherwise, you will have to change the glob pattern in the script.

When all is set, run the script from your new directory

```bash
pwd
$ /Filepath/To/MyNewDirectory
bash ./myrepo/ScerSpliceSeqPipeline.sh
```

## Summary of output

the output directory '/Filepath/To/MyNewDirectory/STAR_AllReadsAlignments' contains alignment files for each pair-end pair of fastq files. These alignment files are for data where PCR duplicates reads are not collapsed.

the  output directory '/Filepath/To/MyNewDirectory/STAR\_DuplicateReadsRemovedAlignments' contains alignment files where reads are collapsed. PCR duplicates Reads are collapsed by filtering out reads (pre-alignment) that are not unique with respect to all base calls. This file also contains a file that  This directory also contains a + strand and a - strand bedgraph that depicts coverage of inserts based on paired-end read locations (note that in calculating this coverage file, there is an underlying assumption that there are no spliced segments in the insert that were not sequenced). The file 'intronCoverageCombined\_SJ_outmerged.tab' file counts number of exon-exon junction (EE), intron-exon junction, and total number of intron-mapping reads for each intron.

The file /Filepath/To/MyNewDirectory/MappingStats.txt is a tab-delimited text file that summarizes for each sample, how many on-target, off-target, and unmappable reads there were for each sample before and after filtering PCR duplicates.

The files in /Filepath/To/MyNewDirectory/Branchpoint_heatplots are plots like below:
![alt text](./images/100K\_Subsampled.BranchPointAbsoluteMeanPeak.pdf)
![alt text](./images/100K\_Subsampled.BranchPointRelativeMeanCutoff.pdf)
![alt text](./images/100K\_Subsampled.TSSAbsoluteMeanPeak.pdf)
![alt text](./images/100K\_Subsampled.TSSRelativeMeanCutoff.pdf)

These plots were made using the SampleData which has only 100K reads. A real dataset would make the plots look better. The underlying matrix of read counts at each position around the branch point for each gene can be found in the '/Filepath/To/MyNewDirectory/STAR\_DuplicateReadsRemovedAlignments' directory. 

