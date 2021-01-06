#!/usr/bin/env nextflow
// @bioinformatic @nextflow @pipeline

// Pipeline to optimise the QC analysis of reads. This is a first trial for
// the analysis of transcriptomics analysis. Here the pipeline focus on the
// session1 of the tuto given by Abhilash Kannan

//                      DEFINE COMMON PARAMETERS and Arguments and variables:

// Define pipeline parameters require the "params.". "params" can be speciedon the command line by prefixingthe parameter name with a double dash character --paranNamesMost of the input variables are included in a config files that must be run with this recipe.

params.Input_Files = '/home/eduranda/Documents/SwitchDrive/CHUV/Sanglard/Data/Sample_dataset/*.fastq.gz'
params.Output_dir = '/home/eduranda/Documents/SwitchDrive/CHUV/Sanglard/Bioinformatics/AK_Tuto_Transcriptomic/20200409_Session1/'

//import variable from the config file
ClumpVar=params.ClumpVar
RmAdaptor_Var=params.Rm_Adaptor
QFilter_Var=params.QFilter
ArtFilter_Var=params.ArtifactFilter
Conta_Removal=params.Conta_Removal
Complex_Filter=params.Complex_Filter
Rm_rRNA=params.Rm_rRNA
Hisat2_var=params.Hisat2
featureCounts_Var=params.featureCounts_Var

// get the species liste
Species=params.Species

// prepare a string that containing the list "-ref <fastafile>" patern that will be included in the bash command.
SortmeRNA_Ref_Pattern="-ref ${Rm_rRNA.Ref_Databases_Name.join(" -ref " )}"

// define the path of the index for the database of
SortmeRNA_DBIndex_Path=Rm_rRNA.DB_Index_Path

// Create variable containing the Output directory
out = params.Output_dir

// get the pattern to get the sample name the Illumina replicates
Pattern_getName_WOilluminaRep=params.Pattern_getSampleName_NoReplicate

// get the list of files
File_List=file(params.Input_Files)

//                        DEFINE CHANNELS :
// Channels are essential to communicate between processes. Here this is a list
// of fastqc files that will to process.

// Channel containing files
File_List_ch = Channel.from(File_List)


// if the Reads has been corrected already fill in the channel used as input for the alignement. Otherwise,
// fill in the channel for starting the correction
SkipCorrection= file("${PWD}/CorrectedReads/rRNA_Removal").exists()

(File_Idx_Init_ch, C_Skip_Corr) = ( SkipCorrection
                 ? [Channel.empty(), Channel
                                       .fromPath("${PWD}/CorrectedReads/rRNA_Removal/No_rRNA*.fastq.gz")
                                       .map { file -> tuple(file.getSimpleName(), file) }]
                 : [Channel
                     .from(File_List)
                     .map { file -> tuple(file.getSimpleName(), file) }, Channel.empty()] )

//Duplicate the channel with the input files liste
File_Idx_Init_ch.into { File_Idx_ch ; File_Idx_QC_ch }

// Prepare channel containing the SortmeRNA database file
Ch_SortmeRNA_DB_File_Address=Channel
  .from(Rm_rRNA.Ref_Databases_Name)
  .map { Rm_rRNA.DB_File_Path + it}


Channel.from(Hisat2_var.GenomeFASTAweb).into { C_Fasta; C_Fasta2}
C_geneID=Channel.from(featureCounts_Var.g)
C_GFF=Channel.from(featureCounts_Var.Web_Features_File)

// prepare channels for the alignement
C_Genomes=Channel
  .from(Hisat2_var.GenomeID)
  .merge(C_Fasta)


C_para_Hisat2=Channel
  .from(Hisat2_var.GenomeID)
  .merge(Channel.from(Hisat2_var.max_intronlen))
  .join(Channel.from(Species))
  .toList()



// C_Genomes.view()
// C_para_Hisat2.view()

C_GFF_n_gID=Channel
  .from(Hisat2_var.GenomeID)
  .merge(C_GFF,C_geneID )
  .join(Channel.from(Species))


C_Sel_Genome=C_Genomes
  .join(Channel.from(Species))


C_Sel_GFF=C_GFF_n_gID
  .join(Channel.from(Species))

// prepare the index folder for the alignement if required
if( file(Hisat2_var.FolderIndx).exists() ){
  println("Folder for Hisat2 Indx exist")
  if (file(Hisat2_var.FolderIndx).list().size()<2){
    newFile=file("$Hisat2_var.FolderIndx/Index.txt")
    newFile.write("Indexing created on ")
  }
  CreatedIndxFolder=Hisat2_var.FolderIndx

}
else{
  println("Creating folder for Hisat2 indexing")
  file(Hisat2_var.FolderIndx).mkdirs()
  CreatedIndxFolder=Hisat2_var.FolderIndx
  newFile=file("$Hisat2_var.FolderIndx/Index.txt")
  newFile.write("Indexing created on ")
}

Channel.fromPath("${CreatedIndxFolder}/*").into {IndxFiles; IndxFilesEmpty }

// build the Hisat2 index if Requires
process Build_index {
  publishDir "$Hisat2_var.FolderIndx" , mode : 'copy' ,pattern : "*.ht2"
  publishDir "$Hisat2_var.FolderIndx" , mode : 'copy' ,pattern : "*.txt"

  input:
  path FIdx from IndxFiles.collect()
  tuple val(x),  path(Fasta) from C_Sel_Genome

  output:
  // file "*"
  path "${x}*.ht2" into C_Index_Ini includeInputs true

  script:
  if((FIdx.name.toString() =~ /${x}/).size()==0)

  // """
  // echo "Gooood ${x} envoie le paté avec ${Fasta}" > output.txt
  // wget $Fasta -O ${x}_Ref.fasta.gz
  // gunzip ${x}_Ref.fasta.gz
  // hisat2-build ${x}_Ref.fasta $x -p 8
  // """
  """
  echo "Gooood ${x} envoie le paté avec ${Fasta}" > output.txt
  gunzip -c $Fasta > ${x}_Ref.fasta
  hisat2-build ${x}_Ref.fasta $x -p 8
  """

  else
  """
  ls . > output.txt
  """
}
// C_Index.view()
C_Index_Ini.into { C_Index; C_Index2}

// Hisat2_var.FolderIndx.view()
process Clupify {
  cache true
  // setup maxmemoryused by the processes
  memory '4 GB'
  // The publishDir allows you to store the process outputs in a directory of your choice.
  // Specify the saveAs parameter to give each file a name of your choice proving a custom
  // rule as a closure
  // !!Produce only a link to the file
  publishDir "$out/Clumpified/"


  // !! produce a hard copy of the file
  // 'publishDir "$PWD/Clumpified_fastQ", mode: "copy" '
  //s
  //import the list of files containing reads to process
  input:
    tuple x , file(Input) from File_Idx_ch
    // val Name from File_names_Clump_ch
  output:
    // file "Started_File.txt"
    tuple x, file("Clumpified_${x}.fastq.gz") into ClumpReads
    tuple x, val(Input) into File_IDX_Name_ch
    // file "Clumpified_${x}.fastq.gz" into ClumpReads

  script:
        // #echo File number $x of name $Name > Started_File.txt
    """
    clumpify.sh -Xmx2g zl=$ClumpVar.zl reorder=$ClumpVar.reorder in=$Input out="Clumpified_${x}.fastq.gz"
    """
  }

process Rm_adapters{
  memory "4 GB"
  publishDir "$out/Adaptor/"

  input:
    tuple x2, file(Input2) from ClumpReads
    // file Input2 from ClumpReads
    // val Name2 from File_names_RmAdap_ch
    path Ref from RmAdaptor_Var.ref
  output:
    tuple x2, file("Rm_Adaptors_${x2}.fastq.gz") into Adapter_Trim
    // file "Rm_Adaptors_${Name2}.fastq.gz" into Adapter_Trim
  script:
    // println { "Trimming adapters from $Name" }
    """
    bbduk.sh -Xmx2g ktrim=$RmAdaptor_Var.ktrim ordered minlen=$RmAdaptor_Var.minlen minlenfraction=$RmAdaptor_Var.minlenfraction mink=$RmAdaptor_Var.mink tbo tpe k=$RmAdaptor_Var.k ow=$RmAdaptor_Var.ow hdist=$RmAdaptor_Var.hdist \
    hdist2=$RmAdaptor_Var.hdist2 ftm=$RmAdaptor_Var.ftm zl=$RmAdaptor_Var.zl in=$Input2 out="Rm_Adaptors_${x2}.fastq.gz" ref=adapters
    """
    // """
    // bbduk.sh -Xmx2g ktrim=$RmAdaptor_Var.ktrim ordered minlen=$RmAdaptor_Var.minlen minlenfraction=$RmAdaptor_Var.minlenfraction mink=$RmAdaptor_Var.mink tbo tpe k=$RmAdaptor_Var.k ow=$RmAdaptor_Var.ow hdist=$RmAdaptor_Var.hdist \
    // hdist2=$RmAdaptor_Var.hdist2 ftm=$RmAdaptor_Var.ftm zl=$RmAdaptor_Var.zl in=$Input2 out="Rm_Adaptors_${Name2}.fastq.gz" ref=$Ref
    // """
}

process QFiltering{
  memory "4 GB"

  publishDir "$out/QFiltering"

  input:
    tuple x3, file(Input3) from Adapter_Trim
    path Ref from QFilter_Var.ref

  output:
    tuple x3, file("QFilter_${x3}.fastq.gz") into QFiltered

  script:
    """
     bbduk.sh -Xmx2g maq=$QFilter_Var.maq trimq=$QFilter_Var.trimq qtrim=$QFilter_Var.qtrim ordered ow=$QFilter_Var.ow maxns=$QFilter_Var.maxns minlen=$QFilter_Var.minlen minlenfraction=$QFilter_Var.minlenfraction \
     k=$QFilter_Var.k hdist=$RmAdaptor_Var.hdist cf=$QFilter_Var.cf ref=artifacts in=$Input3 out="QFilter_${x3}.fastq.gz"
    """
}
// //
process ArtifactFilter{
  memory "4 GB"

  publishDir "$out/Artifact_Filter"

  input:
    tuple x4, file(Input4) from QFiltered
    path Ref from QFilter_Var.ref

  output:
    tuple x4, file("ArtFilter_${x4}.fastq.gz") into WOArtifact_ch

  script:
    """
    bbduk.sh -Xmx2g ordered ow=$ArtFilter_Var.ow k=$ArtFilter_Var.k zl=$ArtFilter_Var.zl hdist=$ArtFilter_Var.hdist in=$Input4 out="ArtFilter_${x4}.fastq.gz" ref=artifacts
    """
}

process Contamiation_Removal {
  memory "4 GB"

  publishDir "$out/Conta_rm"

  input:
    tuple x5, file(Input5) from WOArtifact_ch
    val Ref from Conta_Removal.ref

  output:
    file "Contaminants_stats_${x5}.txt"
    tuple x5, file("NoConta_${x5}.fastq.gz") into WOConta_ch
    file "MatchConta_${x5}.fastq.gz"

  script:
  """
    wget $Ref -O phix174_ill.ref.fa.gz
    bbduk.sh -Xmx2g in=$Input5 out="NoConta_${x5}.fastq.gz" outm=MatchConta_${x5}.fastq.gz ref="phix174_ill.ref.fa.gz" k=$Conta_Removal.k hdist=$Conta_Removal.hdist stats=Contaminants_stats_${x5}.txt
  """
}

process Filter_Low_Complexity {
  memory "4 GB"

  publishDir "$out/Filter_Low_Complexity"

  input:
    tuple x6, file(Input6) from WOConta_ch

  output:
    tuple x6, file("HighComp_${x6}.fastq.gz") into LowComp
    file "LowComp_${x6}.fastq.gz"
  script:
    """
    bbduk.sh in="$Input6" out=HighComp_${x6}.fastq.gz outm=LowComp_${x6}.fastq.gz entropy=$Complex_Filter.entropy
    """
}

process Remove_rRNA{
  memory "12 GB"
  publishDir "$out/rRNA_Removal", mode : 'copy'

  input:
    tuple x7, file(Input7) from LowComp
    path Ref_DB_Index from SortmeRNA_DBIndex_Path
    val Ref_Pattern from SortmeRNA_Ref_Pattern
    // val kvdb from Rm_rRNA.kvdb
    path Ref_Files from Ch_SortmeRNA_DB_File_Address.collect()

  output:
    tuple x7, file("No_rRNA_${x7}.fastq.gz") into NorRNA
    file "With_rRNA_${x7}.fastq.gz"
    file "Rm_rRNA_${x7}.log"
  script:
    """
      mkdir -p kvdb/
      zcat $Input7 > Reads.fastq
      sortmerna $Ref_Pattern -reads Reads.fastq --num_alignments $Rm_rRNA.num_alignments --kvdb kvdb --idx $Ref_DB_Index -threads 4 --fastx --aligned "With_rRNA_${x7}" --other "No_rRNA_${x7}" --paired_in > Rm_rRNA_${x7}.log
      gzip No_rRNA_${x7}.fastq
      gzip With_rRNA_${x7}.fastq
      rm Reads.fastq
      rm -rf kvdb/
    """
}

// duplicate the corrected file channel NorRNA one channel for the QC, one for alignement
NorRNA.into {QC_Corrected_Reads ; Corrected_Reads}

// add content of the "RmRNA" channel to the "File_Idx_QC_ch". Wi do this by combining both a copy of the channel containing all input fastq.gz files with those from the channel obtained from the "Remove_rRNA" channel. we finally group files based on the based on their base on their key (name w/o extension).
QC_Corrected_Reads
  .mix(File_Idx_QC_ch)
  .groupTuple()
  .set { Combined_Reads}

process QC {
  memory "4 GB"
  // publishDir "Results/QC_reports", pattern: '*.html'
  publishDir "Results/QC_reports", pattern: '*.html', mode : 'move'
  publishDir "$out/QC_reports", pattern: '*.zip'
  input:
    tuple x, file(FileList) from Combined_Reads
    // each file(Files) from FileList
  output:
    file "*.zip" into QC_Files_ch
    file "*.html"
  script:
    """
      fastqc $FileList
    """
}

process MultiQC{
  memory "4 GB"
  publishDir "Results/QC_reports", mode: 'move'

  input:
  file Zip_Files from QC_Files_ch.collect()

  output:
  file "*"

  script:
  """
    multiqc $Zip_Files
  """
}

// Group "corrected reads" file by sample. This combine all replicates all together. This is because Hisat2 can combine reads from all fastq files.
Corrected_Reads
  .mix(C_Skip_Corr)
  .map { a, file -> tuple( (file.name.toString() =~ Pattern_getName_WOilluminaRep)[0][1] , file ) }
  .groupTuple()
  .set { Grp_CorrReads }


// Alignement using HiSAT2 -> first species
process Alignement_1{
  // memory "4 GB"
  memory "8 GB"
  cpus 4
  publishDir "$out/Alignement/$Species"


  input:
    // val Species from Hisat2_var.GenomeID[0]
    val Species from C_para_Hisat2.val[0][0]
    val IntroLen from C_para_Hisat2.val[0][1]
    tuple x, file(Files) from Grp_CorrReads
    path ref from C_Index.collect()
    // path Ref_Gen from Hisat2_var.ReferenceGenome_Web

  output:
    // tuple x , file("Alignement_${x}.bam") into Bam_files_ch
    tuple Species,  file("Alignement_${x}.bam") into BamFiles
    tuple x, file("Not_Aligned_1_${x}.fastq.gz") into Align_Again

  script:
    """
      zcat $Files > Test.fastq
      hisat2 -p 4 -x $Species -U Test.fastq --max-intronlen $IntroLen --un-gz Not_Aligned_1_${x}.fastq.gz | samtools sort > Alignement_${x}.bam
    """
}

if (C_para_Hisat2.val.size()>1){
  // Alignement using HiSAT2 -> second
  process Alignement_2{
    memory "8 GB"
    cpus 4
    publishDir "$out/Alignement/$Species"

    input:
      // val Species from Hisat2_var.GenomeID[0]
      val Species from C_para_Hisat2.val[1][0]
      val IntroLen from C_para_Hisat2.val[1][1]
      tuple x, file(Files) from Align_Again
      path ref from C_Index2.collect()
      // path Ref_Gen from Hisat2_var.ReferenceGenome_Web

    output:
      // tuple x , file("Alignement_${x}.bam") into Bam_files_ch
      tuple Species,  file("Alignement_${x}.bam") into BamFiles2
      file "Not_Aligned_1_${x}.fastq.gz" into Align_Again2

    script:
      """
        zcat $Files > Test.fastq
        hisat2 -p 4 -x $Species -U Test.fastq --max-intronlen $IntroLen --un-gz Not_Aligned_1_${x}.fastq.gz | samtools sort > Alignement_${x}.bam
      """
  }
}else{
  Align_Again2=Channel.empty()
  BamFiles2=Channel.empty()
}

Bams=BamFiles
  .mix(BamFiles2)
  .groupTuple()

// C_Sel_GFF.view()

// Counting the Reads mapped to the reference Genome using
process Read_counting{
  memory "10 GB"
  publishDir "Results/ReadCounts/${Species}", mode: 'move' , overwrite: true

  input:
  tuple val(Species), file(Bam) from Bams
  tuple val(SameSpecies), path(GFFs), val(g) from C_Sel_GFF

  output:
  file "*"

  script:
  // """
  // ls > output.txt
  // """
  """
    featureCounts -a $GFFs -F $featureCounts_Var.F -o Read_Count.txt *.bam -t $featureCounts_Var.t -g $g --extraAttributes $featureCounts_Var.extraAttributes \
    -Q $featureCounts_Var.Q > Read_Stats.txt
  """
}
// println(Result)
