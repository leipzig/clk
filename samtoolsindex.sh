set -euo pipefail
cd $TMPDIR

export PATH=/opt/conda/bin:$PATH
/opt/conda/bin/aws s3 cp s3://panorama-clk-repro/${project}/${sample}.Aligned.sortedByCoord.out.bam .

/opt/conda/bin/samtools index ${sample}.Aligned.sortedByCoord.out.bam

/opt/conda/bin/aws s3 cp ${sample}.Aligned.sortedByCoord.out.bam.bai s3://panorama-clk-repro/${project}/