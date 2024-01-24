#5. VCF Polarization
# Download and set up the polarization script
wget https://github.com/kullrich/bio-scripts/blob/master/vcf/polarizeVCFbyOutgroup.py
chmod +x polarizeVCFbyOutgroup.py

# Ensure Python is installed
./polarizeVCFbyOutgroup.py -h
bcftools query -l original_project_data_outgroup_GTonly.vcf.gz

# Polarize the VCF file
# IND -> specify individual idx (an integer giving position of outgroup in the list of samples in the VCF file got by running aobe command) to be used for switch REF and ALT allele
./polarizeVCFbyOutgroup.py -vcf original_project_data_outgroup_GTonly.vcf.gz -ind IND -out original_project_data_outgroup_GTonly_polarized.vcf

# Retain the sample from the original project VCF and remove the outgroup
module unload bcftools
module unload gnu
module load vcftools
vcftools --vcf original_project_data_outgroup_GTonly_polarized.vcf --keep original_project_data_sample_names.txt --recode --stdout | bgzip -c > original_project_data_sample_GTonly_polarized.vcf.gz

# 6. Estimate SFS
# You can now use the resulting file to estimate and plot 1D unfolded SFS and 2D SFS using vcf2sfs.
