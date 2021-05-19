set -euo pipefail
cd $TMPDIR

aws s3 cp s3://clk-splicing/${project}/${sample}_1.fastq.gz .
aws s3 cp s3://clk-splicing/${project}/${sample}_2.fastq.gz .

#https://github.com/ewels/AWS-iGenomes
aws s3 sync s3://clk-splicing/refs/GRCh38_star/ GRCh38_star/


STAR --runMode alignReads \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM ${bytes} \
     --readFilesCommand zcat \
     --outFilterType BySJout   --outFilterMultimapNmax 20 \
     --outFilterMismatchNmax 999   --alignIntronMin 25  \
     --alignIntronMax 1000000   --alignMatesGapMax 1000000 \
     --alignSJoverhangMin 8   --alignSJDBoverhangMin 5 \
     --sjdbGTFfile GRCh38_star/genes.gtf \
     --genomeDir GRCh38_star \
     --runThreadN ${cpus} \
     --outFileNamePrefix ${sample}.  \
     --readFilesIn  ${sample}_1.fastq.gz ${sample}_2.fastq.gz

/code/samtools-1.9/samtools index ${sample}.Aligned.sortedByCoord.out.bam

aws s3 cp ${sample}.Aligned.sortedByCoord.out.bam s3://clk-splicing/${project}/
aws s3 cp ${sample}.Aligned.sortedByCoord.out.bam.bai s3://clk-splicing/${project}/
aws s3 cp ${sample}.Log.final.out s3://clk-splicing/${project}/
aws s3 cp ${sample}.Log.progress.out s3://clk-splicing/${project}/
aws s3 cp ${sample}.SJ.out.tab s3://clk-splicing/${project}/

#STAR defaults
#alignSJDBoverhangMin        3
#alignSJoverhangMin          5
#outFilterMultimapNmax           10
#outFilterMismatchNmax           10
#alignIntronMin              21
#alignIntronMax              0
#alignMatesGapMax            0
#sjdbOverhang                            100

# from lr2rmats
# "{params.star} --runThreadN {threads}  --genomeDir {input.genome} --readFilesIn {input.read1} {input.read2} "
# "--outFileNamePrefix {params.prefix}.STAR --outSAMtype BAM Unsorted "
# "--outFilterType BySJout   --outFilterMultimapNmax 20 "
# "--outFilterMismatchNmax 999   --alignIntronMin 25   --alignIntronMax 1000000   --alignMatesGapMax 1000000 "
# "--alignSJoverhangMin 8   --alignSJDBoverhangMin 5   --sjdbGTFfile {input.gtf}  --sjdbOverhang 100 > {log} 2 >& 1; "

#from shortmap
# "STAR --runThreadN {}  --genomeDir {} --readFilesIn {} --outFileNamePrefix {}/STAR_mapping/samp{}.STAR --outSAMtype BAM SortedByCoordinate " \
#                       "--outFilterType BySJout   --outFilterMultimapNmax 20 --outFilterMismatchNmax 999   --alignIntronMin 25   " \
#                       "--alignIntronMax 1000000   --alignMatesGapMax 1000000 --alignSJoverhangMin 8   --alignSJDBoverhangMin 5  " \
#                       " --sjdbGTFfile {}  --sjdbOverhang 100 > {} 2 >& 1; samtools index {}/STAR_mapping/samp{}.STARAligned.sortedByCoord.out.bam

#from darts
#      --alignEndsType EndToEnd \