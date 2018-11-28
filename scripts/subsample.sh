set -euo pipefail
cd $TMPDIR

export PATH=/opt/conda/bin:$PATH
echo "downloading s3://panorama-clk-repro/${project}/${sample}.Aligned.sortedByCoord.out.bam"
/opt/conda/bin/aws s3 cp s3://panorama-clk-repro/${project}/${sample}.Aligned.sortedByCoord.out.bam .

echo "running /opt/conda/bin/samtools view -s ${subfraction} ${sample}.Aligned.sortedByCoord.out.bam"
/opt/conda/bin/samtools view -s ${subfraction} -b ${sample}.Aligned.sortedByCoord.out.bam  > ${sample}.Aligned.sortedByCoord.out.subsampled.bam

/opt/conda/bin/samtools index ${sample}.Aligned.sortedByCoord.out.subsampled.bam

echo "uploading ${sample}.Aligned.sortedByCoord.out.bam.bai ..."
/opt/conda/bin/aws s3 cp ${sample}.Aligned.sortedByCoord.out.subsampled.bam s3://panorama-clk-repro/${project}/
/opt/conda/bin/aws s3 cp ${sample}.Aligned.sortedByCoord.out.subsampled.bam.bai s3://panorama-clk-repro/${project}/