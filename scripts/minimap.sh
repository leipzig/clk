set -euo pipefail
cd $TMPDIR
aws s3 cp s3://clk-splicing/${project}/${sample}.fastq.gz .
aws s3 sync s3://clk-splicing/refs/GRCh38_minimap/ GRCh38_minimap/

minimap2 -ax splice -ub -t ${cpus} GRCh38_minimap/genome.fa.smmi ${sample}.fastq.gz > ${sample}.sam 2> ${sample}.minimap.log

aws s3 cp ${sample}.sam s3://clk-splicing/${project}/
aws s3 cp ${sample}.minimap.log s3://clk-splicing/${project}/
