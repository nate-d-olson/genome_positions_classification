## processing samtools mpileup vcf files for loading into R

#ran in mpileup directory - note files generated by bpipe script were in bcf not vcf format
for vcf in *vcf; do ../../../../src/samtools-bcftools-htslib-1.0_x64-linux/bin/bcftools view $vcf -O v -o mpileup_vcf/$vcf; done

