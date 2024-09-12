#!/bin/bash
### BEGINNING OF GENOTYPE DECONVOLUTION SCRIPT ###
# AUTHOR: Eddie Cano-Gamez (ecg@sanger.ac.uk)

#source ~/.bashrc
eval "$(conda shell.bash hook)"
export PATH="/software/hgi/installs/anaconda3/condabin:$PATH"
conda activate vireo_cf14_v2

#bamfile
bcftools_path=/lustre/scratch118/humgen/resources/conda_envs/jupyterhub/bin/bcftools 
snplist=${ResourcesFolder}grch38/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz 
genotypes=${ResourcesFolder}/hipsci_genotypes/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.all_chr.hg38.sorted_unique_annotate_MAF5_filterMissing_added_DS.vcf.gz 
cellRangerPath=$OutFolder/1_RunCellRanger/${ID}/ 
bamFile=${cellRangerPath}"possorted_genome_bam.bam"
barcodeFile=${cellRangerPath}"filtered_feature_bc_matrix/barcodes.tsv.gz"

# Defining messages
printHelpMessage() {
	echo ""
	echo ":::::::::::::::::::::::::::: GENOTYPE DECONVOLUTION WITH CELLSNP AND VIREO ::::::::::::::::::::::::::::"
	echo ""
	echo "This script performs deconvolution of single cells by genotypes using the cellSNP and Vireo algorithms"
	echo "This script can either be run on a single sample or parallelized as a job array"
	echo "The script takes as input the name(s) of the sample(s) and a metadata table with donor numbers and donor IDs for each sample."
	echo "For more information, see usage or refer to the README.md file"
	echo ""
	echo "usage: deconvoluteCells.sh [OPTIONAL ARGUMENTS] sampleNames sampleMetadata"
	echo "" 
	echo "Optional arguments:"
	echo "[-h | --help]          Display this help message"
	echo "[-g | --genotypes]     VCF file with reference genotypes, if available (defaults to no reference genotypes)"
	echo "[-t | --genotypeTag]   Format of genotype tags, if reference genotypes are available (i.e. GT, DP, PL; defaults to GT)"
	echo "[-o | --out]			 Output directory where Vireo results should be stored (defaults to the current working directory)"
	echo "[-p | --parallelized]  If this flag is added, the script will run for a group of samples in parallel as a job array"
	echo ""
	echo "Required Arguments:"
	echo "[-s | --snplist]       List of common SNPs used for genotype calling with cellSNP (defaults to a file with hg38 SNPs from 1K genomes)"
  echo "sampleNames         Name of the sample to process (e.g. I1484). If running in --parallelized mode, then a text file with a list of sample names (one per line)"
	echo "sampleMetadata	  Text file containing three columns in the following order: sample IDs, number of donors per sample and sample composition (i.e. donor IDs present in each sample, separated by semicolons)"
	echo ""
	echo "IMPORTANT: Always provide optional arguments before positional arguments"
	echo ""
}

printErrorMessage() {
	echo ""
  	echo "[deconvoluteCells]: ERROR: $1"
  	echo "Please see program usage (--help)"
  	echo ""
}

# Setting help message as default program call
if [ "$#" -eq 0 ]
then
    printHelpMessage
    exit 1
fi


# Setting arguments to default values
genotypeTag="GT"
outputPath="vireo/"
bcftools_path=/lustre/scratch118/humgen/resources/conda_envs/jupyterhub/bin/bcftools
vireo_path=/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/vireo
cellsnp_path=/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/cellSNP_dir/cellSNP-0.1.7/cellSNP/cellSNP_en6.py


#echo $cellRangerPath
# Fetching optional arguments
# TO DO: Add an option to allow keeping intermediary files
for arg in "$@"
do
    case $arg in
        -h|--help)
        printHelpMessage
        exit 1
        ;;
        -s|--snplist)
        snplist="$2"
        shift # Remove --snplist from argument list
        shift # Remove argument value from argument list
        ;;
        -dc|--donor_col)
        donor_col="$2"
        shift # Remove --snplist from argument list
        shift # Remove argument value from argument list
        ;;
        -g|--genotypes)
        genotypes="$2"
        shift # Remove --genotypes from argument list
        shift # Remove argument value from argument list
        ;;
        -t|--genotypeTag)
        genotypeTag="$2"
        shift # Remove --genotypes from argument list
        shift # Remove argument value from argument list
        ;;
        -o|--out)
        outputPath="$2"
        shift # Remove --out from argument list
        shift # Remove argument value from argument list
        ;;
        -ID|--sampleID)
        sampleID="$2"
        shift # Remove --out from argument list
        shift # Remove argument value from argument list
        ;;
        -id|--ID_ind)
        ID_ind="$2"
        shift # Remove --out from argument list
        shift # Remove argument value from argument list
        ;;
        -m|--sampleMetadata)
        sampleMetadata="$2"
        shift # Remove --out from argument list
        shift # Remove argument value from argument list
        ;;
        -c|--cellRangerPath)
        cellRangerPath="$2"
        shift # Remove --out from argument list
        shift # Remove argument value from argument list
        ;;
        -p|--parallelized)
        parallelized="$2"
        shift # Remove --out from argument list
        shift # Remove argument value from argument list
        ;;
        -r|--rerun)
        rerun="$2"
        shift # Remove --out from argument list
        shift # Remove argument value from argument list
        ;;
        -bcf|--bcftools_path)
        bcftools_path="$2"
        shift # Remove --out from argument list
        shift # Remove argument value from argument list
        ;;
    esac
done

# Fetching positional arguments
echo ""
echo "[deconvoluteCells]: Running deconvoluteCells script..."

echo "[deconvoluteCells]: Running script in standard mode..."
  #sampleID="$1"
  echo "[deconvoluteCells]: Processing sample ${sampleID}..."
  

#############################################################################################################
############################################### CELLSNP #####################################################
#############################################################################################################

# Verifying that input files for cellSNP exist
bamFile=${cellRangerPath}"possorted_genome_bam.bam"
barcodeFile=${cellRangerPath}"filtered_feature_bc_matrix/barcodes.tsv.gz"


echo "[deconvoluteCells]: Fetching BAM file..."
echo $bamFile
if [ -f $bamFile ]
then
	:
else
	printErrorMessage "BAM file does not exist [${sampleID}]"
	exit 1
fi

echo "[deconvoluteCells]: Fetching 10X cell barcodes..."
if [ -f $barcodeFile ]
then
	:
else
	printErrorMessage "Barcode file does not exist [${sampleID}]"
	exit 1
fi

echo "[deconvoluteCells]: Fetching SNP list..."
if [ -f $snplist ]
then
	:
else
	printErrorMessage "SNP list file does not exist"
	exit 1
fi


# Creating temporary directory for intermediary files
temporaryDirectory=${outputPath}"/_vireo_temp"

mkdir -p $temporaryDirectory

# Running cellSNP
cellsnpOutput="${temporaryDirectory}/${sampleID}.vcf"

echo "[deconvoluteCells]: Using cellSNP to call genotypes from 10X reads..."
#/software/teamtrynka/conda/trynka-base/bin/cellSNP
echo "bam file: "$bamFile
echo "barcode file "$barcodeFile
echo "snplist: "$snplist
echo "count:"$cellsnpOutput
python3 /lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/cellSNP_dir/cellSNP-0.1.7/cellSNP/cellSNP_en6.py \
    -s $bamFile \
    -b $barcodeFile \
    -o $cellsnpOutput \
    -R $snplist \
    -p 20 \
    --minMAF 0.01 \
   --minCOUNT 20

# Verifying that input files for Vireo exist
cellsnpVCF="${cellsnpOutput}.gz"

echo "[deconvoluteCells]: Fetching cellSNP output..."
if [ -f $cellsnpVCF ]
then
	:
else
	printErrorMessage "cellSNP VCF file does not exist [${sampleID}]"
	echo ""
	echo "[deconvoluteCells]: Stopping now..."
	exit 1
fi

echo "[deconvoluteCells]: Fetching metadata file..."
if [ -f $sampleMetadata ]
then
	:
else
	printErrorMessage "Metadata file does not exist"
	echo ""
	echo "[deconvoluteCells]: Stopping now..."
	exit 1
fi

vireoDirectory=${outputPath}

mkdir -p $vireoDirectory

# Running deconvolution without reference genotypes
echo $genotypes
if [ -z "$genotypes" ]
then
  echo "[deconvoluteCells]: Running deconvolution without reference genotypes..."
  
  # Identifying donors in sample pool 
  echo "[deconvoluteCells]: Fetching number of donors in input sample..."
  numberOfDonors=$(cat "$sampleMetadata" | grep "$sampleID" | awk '{print $7}')
  echo "[deconvoluteCells]: ${numberOfDonors} donors found in sample..."
  
  # Running Vireo
  echo "[deconvoluteCells]: Using Vireo to deconvolute cells... (no genotypes)" 
  /lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/vireo \
   -c $cellsnpVCF \
   -N $numberOfDonors \
   -o $outputDirectory

# Running with reference genotypes  
else
  echo "[deconvoluteCells]: Running deconvolution with reference genotypes..."
  
  echo "[deconvoluteCells]: Fetching reference genotypes..."
  if [ -f $genotypes ]
  then
    :
  else
    printErrorMessage "Reference genotypes file does not exist"
    echo ""
    echo "[deconvoluteCells]: Stopping now..."
    exit 1
  fi
  
  # Identifying donors in sample pool
  detectedDonors="${temporaryDirectory}/${sampleID}_donorNames.txt"
  
  echo "[deconvoluteCells]: Finding donors present in the input sample..."
  cat "$sampleMetadata" | grep "$sampleID," | awk -F ',' '{print $'${donor_col}'}' | tr ';' '\n' > $detectedDonors

  cat $detectedDonors
  
  # Running bcftools
  #echo "[deconvoluteCells]: Using bcftools to subset reference genotypes..."
  
  subsetFile="${temporaryDirectory}/${sampleID}_subset.vcf.gz"
  
  ${bcftools_path} view \
    $genotypes \
    -R $cellsnpVCF \
    -S $detectedDonors \
    -Oz \
    -o $subsetFile
    
  # Running Vireo
  echo "[deconvoluteCells]: Using Vireo to deconvolute cells..."
  
  python3 /lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/vireo_en6 \
    -c $cellsnpVCF \
    -d $subsetFile \
    -t $genotypeTag \
    -o $outputPath/
      
fi


# Removing intermediary files
#echo "[deconvoluteCells]: Cleaning up..."
#rm -r $temporaryDirectory


echo "[deconvoluteCells]: ...DONE!"
echo " "

### END OF GENOTYPE DECONVOLUTION SCRIPT ###







