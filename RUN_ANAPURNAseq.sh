#!/bin/bash

# Choose the path containing raw fastq files
# Path=$(zenity  --file-selection --title="Choose the folder containing raw fastq files" --file-filter=""Downloads" "Desktop"" --directory)
# echo $Path

if ret=`zenity  --file-selection --title="Choose the folder containing raw fastq files" \
  --filename="/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/Fastq_gz/"  --directory`
  then
    Path=$ret
  else
    zenity --error --text "You must choose the folder where fastq files are located. Start again." --width=300 --height=300
    exit
  fi



#Choose genomes, you can choose until two genomes from the list :["Ca", "Cg", "Mm","Hs", "Mc", "Rm"]

if ret=`zenity --list --width=1000 --height=300 \
  --radiolist --title="Choose your genome of interest" \
  --column="Choice" --column="ID" --column="Specie" --column="Strain" --column="Genome File Name" --column="GFF file Name" \
    "" "Ca" "Candida albicans" "SC5314" "C_albicans_SC5314_A22_current_chromosomes_A.fasta.gz" "http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_A22_current_features.gff" \
    "" "Cg" "Candida glabrata" "DSY562" "DSY562_corrected_genome_03042017_J.fasta.gz" "DSY562_corrected_genome_03042017_J.gff" \
    "" "Mm" "Mus musculus" "GRCm38.p6" "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz" \
    "" "Hs" "Homo sapiens" "GRCh38" "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz" "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz" \
    "" "Mc" "Mucor circinelloides" "1006Phmd" "Mucor_circinelloides_1006Phmod.fasta" "Mucor_circinelloides_1006Phmod.gff" \
    "" "Rm" "Rhizopus microscopus" "ATCC52813" "Rhizopus_microsporus_ATCC52813mod.fasta" "Rhizopus_microsporus_ATCC52813mod.gff"`

    then
    Specie1=$ret
    Specie2=$(zenity --list --width=1000 --height=300\
      --radiolist --title="Choose your genome of interest" \
      --column="Choice" --column="ID" --column="Specie" --column="Strain" --column="Genome File Name" --column="GFF file Name" \
        "" "" "None" "" "" "" \
        "" "Ca" "Candida albicans" "SC5314" "C_albicans_SC5314_A22_current_chromosomes_A.fasta.gz" "http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_A22_current_features.gff" \
        "" "Cg" "Candida glabrata" "DSY562" "DSY562_corrected_genome_03042017_J.fasta.gz" "DSY562_corrected_genome_03042017_J.gff" \
        "" "Mm" "Mus musculus" "GRCm38.p6" "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz" \
        "" "Hs" "Homo sapiens" "GRCh38" "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz" "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz" \
        "" "Mc" "Mucor circinelloides" "1006Phmd" "Mucor_circinelloides_1006Phmod.fasta" "Mucor_circinelloides_1006Phmod.gff" \
        "" "Rm" "Rhizopus microscopus" "ATCC52813" "Rhizopus_microsporus_ATCC52813mod.fasta" "Rhizopus_microsporus_ATCC52813mod.gff")
        cd $Path

        if test -f *.fastq
        then
          zenity --error --text "The folder does not contains fastq files. Start again." --width=300 --height=300
          exit
        else
          #NextFlow File
          Nf_File="/media/edurandau/DATA/Eric_Durandau/Sanglard/Labbook/Workbench/Bioinformatics/ANAPURRNA-seq/20200420-Nextflow_ReadsCorrection_v0.nf"

          #NextFlow config file
          Cf_File="/media/edurandau/DATA/Eric_Durandau/Sanglard/Labbook/Workbench/Bioinformatics/ANAPURRNA-seq/ANAPURNA_Config_.config"

          nextflow $Nf_File -c $Cf_File --Specie1 $Specie1 --Specie2 $Specie2 -resume
        fi

      else
        zenity --error --text "You must choose the first genome for the alignment. Start again." --width=300 --height=300
        exit
      fi
