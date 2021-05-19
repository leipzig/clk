set -euo pipefail
cd $TMPDIR
export PATH=/opt/conda/bin:$PATH

/opt/conda/bin/aws s3 cp s3://clk-splicing/${project}/${sample}.bam .
/opt/conda/bin/aws s3 cp s3://clk-splicing/${project}/${sample}.bam.bai .
mkdir GRCh38_star
/opt/conda/bin/aws s3 cp s3://clk-splicing/refs/GRCh38_star/${gtf} GRCh38_star/


/opt/conda/bin/IsoModule GRCh38_star/${gtf} ${sample}.bam -i 3 -T 10 -c 1 -v 3 -e 50 -C 20 -o .

/opt/conda/bin/aws s3 cp ${sample}.bam.IsoMatrix s3://clk-splicing/${project}/
/opt/conda/bin/aws s3 cp ${sample}.bam.IsoExon s3://clk-splicing/${project}/

