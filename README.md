# ANAPURNA-seq

## 0. Quick start

ANAPURNA-seq is designed to run on a UNIX-based system and **requires [miniconda][5], [zenity][4], Bash 3.2 (or later) and [Java 8 (or later, up to 15)][6]**. The computer (or sever) running ANAPURNA-seq must also be connected to the internet.

Start ANAPURNA-seq by running the ` $ sh RUN_ANAPURNAseq.sh`. A sequence of windows will appear prompting the user to choose reference genomes.

## 1. Introduction

![ANAPURNA](Icons/ANAPURNA-seq.png)

ANAPURNAseq (Automated Nextflow Alignment Pipeline for Unprocessed RNA-seq) perform first, various read correction step, then align corrected reads to one or two reference genomes and finally count the number of reads that can be associated to a chosen feature (often *gene* or *CDS*). The result of a run is basically one or multiple "readcount.txt" files (see *[fig A.1](#fA.1)*).	

<a name="fA.1"></a>![read-procession-steps](Icons/20201209-Scheme-ANAPURNA.png)

*fig A.1 : The various steps performed by ANAPURNA-seq* 

ANAPURNA-seq is using various  software from the Bioconda and Conda Forge channels. All softwares are listed in [Table 1.1](#table1.1) with their corresponding function into the pipeline. The first step consists in compressing fastq.gz files by reordering them by similarity. This is achieved using the clumpify function from the [*bbmap*][9][^9] packages. This has the advantage to improve the speed of the downstream steps. Following compression, extremities of reads are trimmed to remove adaptors. 

Multiple steps of filtering are then applied on the compressed reads: Quality filtering, contamination removal, low complexity filtering, removal of artefacts. These are achieved using again the [*bbmap*][9][^9] packages. The last filtering step is the removal of rRNA using the *[sortmerna][10]*[^10] package .

After reads have been compressed, filtered and corrected, they are ready to be aligned to one or two reference genomes. This is achieved by using first the HiSAT [11][12] software  and second the "featurecount" function of the [*subread*][12] [^12]package . The current list of available genomes are shown in [Table 1.2](#table1.2). At the end, of the whole procedure, QC reports and readcount table(s) are generated. QC report enables to check the quality of the reads pre- and post-processing. Readcount tables can then be used for different downstream analysis such as differential gene expression analysis. 

All the steps have been embedded in a *Nextflow* pipeline. Nextflow enables scalable and reproducible scientific workflows using software containers. It allows the adaptation of pipelines written in the most common scripting languages. Nextflow is based on the *data flow* programming model which greatly simplifies writing complex distributed pipelines. Parallelisation is implicitly defined by the processes input and output declarations. The resulting applications are inherently parallel and can scale-up or scale-out, transparently, without having to adapt to a specific platform architecture. This make the language perfectly suited to run on cluster architectures.

## 2. Procedure to run ANAPURNA-seq 

1. Start ANAPURNA-seq by running the ` $ sh RUN_ANAPURNAseq.sh`
2. Go to the folder where *fastq.gz* have been saved
3. Select the first genome on which reads have to be aligned
4. Select the second genome on which reads, not aligned to the first genome, have to be aligned. Following that step, a terminal window should appear and give you the progress of the processing.
5. When the analysis is done the new message should appear prompting the user to remove or not temporary files. Click OK if you are sure everything has run properly otherwise, re-run ANAPURNA with the same setting.
6. A last window should open indicating in which folder Results can be found

## 3. Notes

* ANAPURNA-seq can be run stop at any time by selecting the terminal window showing the run progress and by clicking Ctrl+C
* As long as the "work" directory is not deleted from the working directory, all processes already completed are saved. This is useful in case the user wish to run the pipeline multiple times on the same dataset. For example, if the users want to run ANAPURNA a second time and align reads on another genome, none of the steps prior alignment will be performed again. Only the alignment and the featurecount steps. This leads to a large gain of time for the other runs.



## Table

### Table 1.1 : List of packages used by ANAPURNA-seq  <a name="table1.1"></a>

| Name      | Version                    | Function                                                     | Reference     |
| --------- | -------------------------- | ------------------------------------------------------------ | ------------- |
| Nextflow  | 21.05.1.5556               | Parallelisation of processes                                 | [3][3][^3]    |
| Zenity    | 3.32.0                     | GUI                                                          | [4][4][^4]    |
| miniconda | 5.10.1                     | Package management system and environment management system  | [5][5][^5]    |
| Java      | openjdk 11.0.11 2021-04-20 | Development environment for building applications and components | [6][6][^6]    |
| fastqc    | 0.11.9                     | reads QC report                                              | [7][7][^7]    |
| multiqc   | 1.10.1                     | combining multiple QC reports in a single file               | [8][8][^8]    |
| bbmap     | 38.90                      | Removal of adaptors and contamination. Quality, artefacts, and low complexity removal. | [9][9][^9]    |
| sortmerna | 5.3.3                      | Removal of rRNA from various species                         | [10][10][^10] |
| samtools  | 1.12                       | Conversion of files from alignment in *.bam*                 | [13][13][^13] |
| hisat2    | 2.2.1                      | Alignment of reads to a genome                               | [11][11][^11] |
| subread   | 2.0.1                      | Counting the number of reads associated to a feature (*featurecount*) | [12][12][^12] |

### Table 1.2: List of genomes currently available in ANAPURNA-seq <a name="table1.2"></a>

| Specie               | Strain               |
| -------------------- | -------------------- |
| *Candida albicans*   | SC5314 [14][14][^14] |
| Candida glabrata     | DSY562               |
| Mus musculus         | GRCm38 [15][15][^15] |
| Homo sapiens         | GRCh38 [16][16][^16] |
| Mucor circinelloides | 1006Phmd             |
| Rhizopus microsporus | ATCC52813            |
| Candida lusitaniae   | DSY4606_V37          |

## 4. References

[1]: https://github.com/Navrique/ANAPURNAseq	"ANAPURNAseq GitHUB"
[2]: ANAPURNA-seq	"ANAPURNA-seq hard copy"
[3]: https://www.nextflow.io/	"nextflow"
[4]: https://doc.ubuntu-fr.org/zenity	"Zenity"
[5]: https://docs.conda.io/en/latest/miniconda.html	"miniconda"
[6]: http://www.oracle.com/technetwork/java/javase/downloads/index.html
[7]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/	"fastqc"
[8]: https://multiqc.info/	"multiqc"
[9]: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/	"bbmap"
[10]: https://github.com/biocore/sortmerna	"sortmerna"
[11]: http://daehwankimlab.github.io/hisat2/	"hisat2"
[12]: http://subread.sourceforge.net/ "subread"
[13]: http://www.htslib.org/	"samtools"
[14]: http://www.candidagenome.org	"Candida genome database"
[15]: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz	"Mouse reference genome"
[16]: ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz	"Human reference genome"
[17]: https://kasperdanielhansen.github.io/genbioconductor/html/limma.html	"limma"

[^1]: https://github.com/Navrique/ANAPURNAseq
[^2]: ./ANAPURNA-seq
[^3]: https://www.nextflow.io/
[^4]: https://doc.ubuntu-fr.org/zenity
[^5]: https://docs.conda.io/en/latest/miniconda.html
[^6]:http://www.oracle.com/technetwork/java/javase/downloads/index.html
[^7]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[^8]: https://multiqc.info/ 
[^9]: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
[^10]: https://github.com/biocore/sortmerna
[^11]: http://daehwankimlab.github.io/hisat2/
[^12]: http://subread.sourceforge.net/
[^13]: http://www.htslib.org/ 
[^14]: http://www.candidagenome.org 
[^15]: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
[^16]: ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz 
[^17]: https://kasperdanielhansen.github.io/genbioconductor/html/limma.html 
[^18]: log2 fold change
[^19]: https://en.wikipedia.org/wiki/Canberra_distance

[20]: http://geneontology.org/	"GO"

[^GO]: http://geneontology.org/
[^string]: https://string-db.org/

[string]: https://string-db.org/	"STRINGdb"
[CP]: https://guangchuangyu.github.io/software/clusterProfiler/	"ClusterProfiler"

[^CP]: https://guangchuangyu.github.io/software/clusterProfiler/



[AnnotationForge]: https://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/SQLForge.pdf	"Annotation Forge"
[NCBI]: https://www.ncbi.nlm.nih.gov/gene/?term=txid5478Organism:exp	"Gene list NCBI"
[ChrFeatures]: http://www.candidagenome.org/download/chromosomal_feature_files/C_glabrata_CBS138/C_glabrata_CBS138_current_chromosomal_feature.tab	"Chromosomal features"
[GeneAssociation]: http://www.candidagenome.org/download/go/	"GO association"
[stringR]: https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html	"STRINGdb package"

[^stringR]: https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html 
[^GeneAssociation]: http://www.candidagenome.org/download/go/
[^ChrFeatures]: http://www.candidagenome.org/download/chromosomal_feature_files/C_glabrata_CBS138/C_glabrata_CBS138_current_chromosomal_feature.tab
[^NCBI]:https://www.ncbi.nlm.nih.gov/gene/?term=txid5478
[^AnnotationForge]: https://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/SQLForge.pdf