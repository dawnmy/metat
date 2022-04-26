
# Kraken standard database ############################################

KRAKEN=/vol/projects/dzhiluo/tools/kraken_uv2000/kraken
KRAKEN_R=/vol/projects/dzhiluo/tools/kraken_uv2000/kraken-report
KRAKEN_MR=/vol/projects/dzhiluo/tools/kraken_uv2000/kraken-mpa-report
db=${1%/}
outdir=${2%/}
NCORE=100
for file in *_cl_mRNA_R1.fastq
do
        sample=${file%_cl_mRNA_R1.fastq}
        outpath=${outdir}/${sample}
        r1=${file}
        r2=${sample}_cl_mRNA_R2.fastq
        logfile=${outpath}_kraken.log
        $KRAKEN --db $db --threads $NCORE --fastq-input --paired $r1 $r2 --unclassified-out ${outpath}_unclassified.fastq --classified-out ${outpath}_classified.fastq --output ${logfile}
        $KRAKEN_R --db $db $logfile > ${logfile%.log}_normal.report
        $KRAKEN_MR --show-zeros --db $db $logfile > ${logfile%.log}_report_mpa.txt
done

# BWA to gene catalog #################################################
export BWA=/vol/biotools/bin/bwa
export NCORE=20
export SAMTOOLS=/vol/biotools/bin/samtools
index=$1
outdir=${2%/}
for file in *_cl_mRNA_R1.fastq
do
        sample=${file%_cl_mRNA_R1.fastq}
        outpath=${outdir}/${sample}
        r1=${file}
        r2=${sample}_cl_mRNA_R2.fastq
        bamfile=${outpath}_sorted.bam
        logfile=${outpath}_aln.log
        samcountfile=${outpath}_samcount.txt
        countfile=${outpath}_count.txt
        $BWA mem -k 31 -t $NCORE $index ${r1} ${r2} | $SAMTOOLS view -Shb - 2>> $logfile|$SAMTOOLS sort -@ 20 -m 10G -f - $bamfile >> $logfile 2>&1
        $SAMTOOLS index $bamfile
        $SAMTOOLS idxstats $bamfile > $samcountfile
        awk -F"\t" '{print $3}' $samcountfile > $countfile

done


