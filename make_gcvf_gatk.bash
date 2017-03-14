fastq="/home/owens/working/Ana/fastq"
ngm="/home/owens/bin/NextGenMap-0.4.12/bin/ngm-0.4.12/ngm"
ref="/home/owens/ref/HA412.v1.1.bronze.20141015.fasta"
bam="/home/owens/working/Ana/bam"
picardtools='/home/owens/bin/picard-tools-2.1.0'
log="/home/owens/working/Ana/log"
project="Ana"
home="/home/owens/working/Ana"
bin="/home/owens/bin"
javarules="-Djava.io.tmpdir=/media/owens/speedy/tmp"
demultiplex='/home/owens/bin/GBS_fastq_Demultiplexer_v8.GO.pl'
plate="Ana"
trim="/home/owens/bin/Trimmomatic-0.32"
trimmeddata="/home/owens/working/Ana/trimmed_data"
unpaired="/home/owens/working/Ana/trimmed_data_unpaired"
ncores="5"
Gpath='/home/owens/bin'
gvcf='/home/owens/working/Ana/gvcf/'
#Make bam.list for GATK
ls -d $bam/*.* | grep bam  | grep -v bai> $home/bamlist.$project.list

#identify local indels

java -Xmx4g -jar $Gpath/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R $ref \
   -I $home/bamlist.$project.list \
   -nt $ncores \
   -log $log/$project.RealignerTargetCreator.log \
   -o $home/$project.realign.intervals

#Realign around local indels
while read prefix
do
java -Xmx4g -jar $Gpath/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $ref \
	-I $bam/$prefix.bam \
	-targetIntervals $home/$project.realign.intervals \
	-o $bam/$prefix.realign.bam \
	-log $log/$plate.$prefix.IndelRealigner.log 

#Call GATK HaplotypeCaller
java -Xmx18g -jar $Gpath/GenomeAnalysisTK.jar \
	-l INFO \
	-R $ref \
	-log $log/$plate.$prefix.HaplotypeCaller.log \
	-T HaplotypeCaller \
	-I  $bam/$prefix.realign.bam \
	--emitRefConfidence GVCF \
	--max_alternate_alleles 2 \
	-variant_index_type LINEAR \
	-variant_index_parameter 128000 \
	-o $gvcf/${prefix}.GATK.gvcf.vcf
done < $home/Samplelist.${project}.txt
exit
