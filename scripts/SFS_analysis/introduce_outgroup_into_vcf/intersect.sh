### Intersect and Merge VCF Files
module load bcftools

#make the names in the vcf-files the same

$ gunzip -c capsella_rubella.vcf.gz > temp.vcf $ awk -F'\t' 'BEGIN {OFS="\t"} {if ($1 ~ /^##/) {print} else if ($1 ~
/^#CHROM/) {gsub("NC_000932.1", "Pltd", $0); gsub("NC_003070.9", "1", $0); gsub("NC_003071.7", "2", $0); gsub("NC_003074.8", "3", $0); gsub("NC_003075.7", "4", $0); gsub("NC_003076.8", "5", $0); gsub("NC_037304.1", "MT", $0); print} else if ($1 !~ /^#/) {gsub("NC_000932.1", "Pltd", $1); gsub("NC_003070.9", "1", $1); gsub("NC_003071.7", "2", $1); gsub("NC_003074.8", "3", $1); gsub("NC_003075.7", "4", $1); gsub("NC_003076.8", "5", $1); gsub("NC_037304.1", "MT", $1); print} else {print}}' temp.vcf > capsella_rubella_renamed.vcf

$ cat capsella_rubella_renamed.vcf | awk '!/^#/ {print $1}' | sort -T /scratch/tali/tmp | uniq
1
2
3
4
5
MT
Pltd

$ bgzip -c capsella_rubella_renamed.vcf > capsella_rubella_renamed.vcf.gz $ bcftools index capsella_rubella_renamed.vcf.gz $ bcftools index arabidopsis_thaliana_whole.genome.vcf.gz
$ bcftools isec -n=2 -p
/scratch/QBIO/student_project_data/Group_4/sfs/bcftools_isec/
arabidopsis_thaliana_whole.genome.vcf.gz capsella_rubella_renamed.vcf.gz $ grep -v "^#" bcftools_isec/0000.vcf | wc -l
112924
$ grep -v "^#" bcftools_isec/0001.vcf | wc -l
112924

#now do the actual intersection

bcftools index arabidopsis_thaliana_whole.genome.vcf.gz
bcftools index capsella_rubella.vcf.gz
bcftools isec -n=2 -p isec arabidopsis_thaliana_whole.genome.vcf.gz capsella_rubella.vcf.gz

# This will generate two files in the 'isec': "0000.vcf" and "0001.vcf"

# Check the number of variant sites in both files, should get the same number for both
grep -v "^#" isec/0000.vcf | wc -l
grep -v "^#" isec/0001.vcf | wc -l

# Identify the files: original VCF with project samples and VCF with outgroup
bcftools query -l isec/0000.vcf
bcftools query -l isec/0001.vcf

# Compress and rename both files

#instead of bgzip you can use    bcftools view file.vcf -Oz -o file.vcf.gz

bgzip -c isec/0000.vcf > original_project_data_isec.vcf.gz
bgzip -c isec/0001.vcf > outgroup_isec.vcf.gz


# Index both compressed files
bcftools index original_project_data_isec.vcf.gz
bcftools index outgroup_isec.vcf.gz

# Merge both files
bcftools merge original_project_data_isec.vcf.gz outgroup_isec.vcf.gz > original_project_data_outgroup.vcf
bgzip -c original_project_data_outgroup.vcf > original_project_data_outgroup.vcf.gz

# Validate the number of sites
zgrep -v "^#" original_project_data_outgroup.vcf.gz | wc -l

# Annotate and extract GT information
bcftools annotate -x ^INFO/GT,^FORMAT/GT original_project_data_outgroup.vcf.gz | bcftools +setGT -- -ta -nu > original_project_data_outgroup_GTonly.vcf
bgzip -c original_project_data_outgroup_GTonly.vcf > original_project_data_outgroup_GTonly.vcf.gz
zgrep -v "^#" original_project_data_outgroup_GTonly.vcf.gz | wc -l