##################################################
## 4. EN‐TEx ATAC‐seq data: downstream analyses ##
##################################################




# metadata file.
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_slims=DNA+accessibility&assay_title=ATAC-seq&biosample_ontology.term_name=stomach&biosample_ontology.term_name=sigmoid+colon&type=Experiment"

# Download peak calling ##

## Create the same folders as ChIP-seq.
mkdir analyses annotation data data/bigBed.files data/bed.files

## save a bigBed with the peaks of interest.
cd analyses
grep -F ATAC-seq ../metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > bigBed.peaks.ids.txt

## Download the file.
cut -f1 bigBed.peaks.ids.txt |\
while read filename; do
  wget -P ../data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

## Check the integrity of our bigBed files with md5sum.
../bin/selectRows.sh <(cut -f1 analyses/bigBed.peaks.ids.txt) metadata.tsv | cut -f1,46 > data/bigBed.files/md5sum.txt

#Go to data and check
cat bigBed.files/md5sum.txt |\
while read filename original_md5sum; do 
  md5sum bigBed.files/"$filename".bigBed |\
  awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
done > tmp 
mv tmp bigBed.files/md5sum.txt

awk '$2!=$3' data/bigBed.files/md5sum.txt
## No output, then it is correct.


#Let's go to annotation directory.
wget -P annotation "https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz"

## Let's uncopress the gtf.gz file.
gunzip gencode.v24.primary_assembly.annotation.gtf.gz


## Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. 

awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed 

## Get from the repository all the TSS and put it no a file. 
wget -P annotation/ "https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/gencode.v24.protein.coding.non.redundant.TSS.bed"


## Convert .bigBed into .bedfiles.
cd analyses
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done

## Intersect the peaks and no-redundant TSS.
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -wa -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed |\
  sort -u > analyses/promoters.peaks."$tissue".ATAC-seq.bed
done


## Overlapping peak.
cut -f2 analyses/bigBed.peaks.ids.txt |\
while read tissue; do
  wc -l analyses/promoters.peaks."$tissue".ATAC-seq.bed
done
# RESULTS:
# 47871 promoters.peaks.sigmoid_colon.ATAC-seq.bed
# 44749 promoters.peaks.stomach.ATAC-seq.bed
# 92620 total


## Intersection of the outer peaks. 
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -v -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.gene.body.bed |\
  sort -u > analyses/non.gene.peaks."$tissue".ATAC-seq.bed
done

## Overlapping peak.
cut -f2 analyses/bigBed.peaks.ids.txt |\
while read tissue; do
  wc -l analyses/non.gene.peaks."$tissue".ATAC-seq.bed
done
# RESULTS:
# 37035 non.gene.peaks.sigmoid_colon.ATAC-seq.bed
# 34537 non.gene.peaks.stomach.ATAC-seq.bed
# 71572 total
		

##################################
## 5.Distal regulatory activity ##
##################################

#Task 1:Create regulatory_elements folder.
mkdir regulatory_elements
