combinedGVCFs='/home/owens/working/Ana/gvcf'
ls $combinedGVCFs | grep "vcf" | grep -v ".idx"   > GVCFs.samplelist.txt
tmp=""
while read prefix
do
        tmp="$tmp --variant $combinedGVCFs/$prefix"
done < GVCFs.samplelist.txt

java -Xmx20g -Djava.io.tmpdir=/home/owens/SB/tmp -jar /home/owens/bin/GenomeAnalysisTK.jar \
        -nt 8 \
        -l INFO \
        -R /home/owens/ref/HA412.v1.1.bronze.20141015.fasta \
        -log log/GenotypeGVCFs.log \
        -T GenotypeGVCFs \
        $tmp \
        -o Ana.GATK.2016.plusUSA.vcf \
	-hets 0.01 \
        --max_alternate_alleles 4
