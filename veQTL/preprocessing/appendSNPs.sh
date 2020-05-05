## This loop makes a list of commands to run and saves to file Append.sh
## It takes the VCF files which have been split by chr they are named like --> largevcf.chr1.recode.vcf only chr changes
## Using those file names I create the file names of the transposed 012 matrcies. These were created by the Matrix012Loop.sh and transposed with Transpose.R. File names --> phg000830_635Indi_chr10.012Mat.txt
## The outfile will contain 5 columns of snp identifies from the vcf file and the 012 matrix for each samples (by column). This file will be saved in the directory it is run, with the same name as the original 012 matrix. 
## WILL OVERWRITE ORIGINAL TRANSPOSED MATRIX 
for f in /STORAGE/George/GTEx/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/split_by_chr/*.vcf
do  
 other="${f/largevcf./Matrix012/phg000830_635Indi_}"
 other2="${other/recode*/012Mat.txt}"
 out="${other/recode*/012Mat.txt}"
 out="${out/*Matrix012\//}"
 printf "grep -v '##' $f | cut -f-5 >tmp.txt\npaste tmp.txt $other2 > $out\n" >> Append.sh
done

## Make the file executable
chmod 755 Append.sh
## runs list of commands (comment out if testing loop etc)
./Append.sh