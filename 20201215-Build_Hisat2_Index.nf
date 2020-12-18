

IndexFolder=params.FolderIndx

Fasta=Channel.from(params.GenomeFASTAweb)
GFF=Channel.from(params.GenomeGFFweb)

C_Genomes= Channel
  .from(params.GenomeID)
  .merge(Fasta)


C_GFF=Channel
  .from(params.GenomeID)
  .merge(GFF)


ToIdx=["Cg", "Ca"]
// ToIdx=["Ca"]
C_Sel_Genome=C_Genomes
  .join(Channel.from(ToIdx))


C_Sel_GFF=C_GFF
  .join(Channel.from(ToIdx))

if( file(IndexFolder).exists() ){
  println("Folder for Hisat2 Indx exist")
  if (file(IndexFolder).list().size()<2){
    newFile=file("$IndexFolder/Index.txt")
    newFile.write("Indexing created on ")
  }
  CreatedIndxFolder=IndexFolder

}
else{
  println("Creating folder for Hisat2 indexing")
  file(IndexFolder).mkdirs()
  CreatedIndxFolder=IndexFolder
  newFile=file("$IndexFolder/Index.txt")
  newFile.write("Indexing created on ")
}

Channel.fromPath("${CreatedIndxFolder}/*").into {IndxFiles; IndxFilesEmpty }


process Build_index {
  publishDir "$IndexFolder" , mode : 'copy' ,pattern : "*.ht2"
  publishDir "$IndexFolder" , mode : 'copy' ,pattern : "*.txt"

  input:
  path FIdx from IndxFiles.collect()
  tuple val(x),  val(Fasta) from C_Sel_Genome

  output:
  // file "*"
  path "${x}*.ht2" into C_Index includeInputs true

  script:
  if((FIdx.name.toString() =~ /${x}/).size()==0)

  """
  echo "Gooood ${x} envoie le paté avec ${Fasta}" > output.txt
  wget $Fasta -O ${x}_Ref.fasta.gz
  gunzip ${x}_Ref.fasta.gz
  hisat2-build ${x}_Ref.fasta.gz $x
  """

  else
  """
  ls . > output.txt
  """
}
C_Index.view()
