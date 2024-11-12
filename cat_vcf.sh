
interactive -p nocona

# make a single vcf for PCA,ADMIXTURE

grep "#" scaffold0001.recode.vcf > pca_20kbp.all.vcf

for i in $( ls *recode.vcf ); do
grep -v "#" $i >> pca_all.vcf;
done

rm *recode*

rm *idx
