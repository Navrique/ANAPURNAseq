#!/bin/bash

# Choose the path containing raw fastq files
# Path=$(zenity  --file-selection --title="Choose the folder containing raw fastq files" --file-filter=""Downloads" "Desktop"" --directory)
# echo $Path
# source ~/.bashrc
echo $PATH

if ret=`zenity  --file-selection --title="Choose the folder containing raw fastq files" \
  --filename="~"  --directory`
  then
    Path=$ret
  else
    zenity --error --text "You must choose the folder where fastq files are located. Start again." --width=300 --height=300
    exit
  fi

# Get the base directory
BASEDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#Choose genomes, you can choose until two genomes from the list :["Ca", "Cg", "Mm","Hs", "Mc", "Rm"]

if ret=`zenity --list --width=1000 --height=300 \
  --radiolist --title="Choose your genome of interest" \
  --column="Choice" --column="ID" --column="Specie" --column="Strain" --column="Genome File Name" --column="GFF file Name" \
    "" "Ca" "Candida albicans" "SC5314" "C_albicans_SC5314_A22_current_chromosomes_A.fasta.gz" "http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_A22_current_features.gff" \
    "" "Cg" "Candida glabrata" "DSY562" "DSY562_corrected_genome_03042017_J.fasta.gz" "DSY562_corrected_genome_03042017_J.gff" \
    "" "Mm" "Mus musculus" "GRCm38.p6" "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz" \
    "" "Hs" "Homo sapiens" "GRCh38" "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz" "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz" \
    "" "Mc" "Mucor circinelloides" "1006Phmd" "Mucor_circinelloides_1006Phmod.fasta" "Mucor_circinelloides_1006Phmod.gff" \
    "" "Rm" "Rhizopus microscopus" "ATCC52813" "Rhizopus_microsporus_ATCC52813mod.fasta" "Rhizopus_microsporus_ATCC52813mod.gff" \
    "" "Af" "Aspergillus fumigatus" "Af_293B" "Af_293B.fasta" "Af_293B.gff" \
    "" "Cl" "Candida lusitaniae" "DSY4606_V37" "Candida_lusitaniae_DSY4606_V37.fasta.gz" "Candida_lusitaniae_DSY4606_V37.gff.gz"`

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
          "" "Rm" "Rhizopus microscopus" "ATCC52813" "Rhizopus_microsporus_ATCC52813mod.fasta" "Rhizopus_microsporus_ATCC52813mod.gff" \
          "" "Af" "Aspergillus fumigatus" "Af_293B" "Af_293B.fasta" "Af_293B.gff" \
          "" "Cl" "Candida lusitaniae" "DSY4606_V37" "Candida_lusitaniae_DSY4606_V37.fasta.gz" "Candida_lusitaniae_DSY4606_V37.gff.gz")
          cd $Path

          FastqGz=$(find . -type f -name "*.fastq.gz")
          # if test -f *.fastq.gz ; then
          if ! [[ (-n $FastqGz)]] ; then
            zenity --error --text "The folder does not contains fastq files. Start again." --width=300 --height=300
            exit
          else
            echo $BASEDIR
            #NextFlow File
            Nf_File="$BASEDIR/20200420-Nextflow_ReadsCorrection_v0.nf"

            #NextFlow config file
            Cf_File="$BASEDIR/ANAPURNA_Config_.config"

            # Start ANAPURNAseq
            $BASEDIR/nextflow $Nf_File -c $Cf_File --Specie1 $Specie1 --Specie2 $Specie2 -resume

            # Message to suppress the work directory
            if zenity --question --text="ANAPURNAseq have ran successfully. Do you want to supress the 'work' directory to save space?"  --width=300 --height=300
              then
                  rm -rf work
                  rm -rf .nextflow*
              fi
            zenity --info --text="Find the result in $Path"  --width=300 --height=300
          fi
    else
      zenity --error --text "You must choose the first genome for the alignment. Start again." --width=300 --height=300
      exit
    fi
