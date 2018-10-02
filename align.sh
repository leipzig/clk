set -euo pipefail
cd $TMPDIR
aws s3 cp s3://{bucket}/{sample}_r1.fq .
aws s3 cp s3://{bucket}/{sample}_r2.fq .
aws s3 sync s3://{bucket}/assets/{reference} .
bwa mem -t {cpus} {reference}.fa {sample}_r1.fq {sample}_r2.fq \
      | samtools sort -o {sample}.bam
samtools index {sample}.bam
aws s3 cp {sample}.bam s3://{bucket}/
aws s3 cp {sample}.bam.bai s3://{bucket}/