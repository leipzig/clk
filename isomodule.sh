set -euo pipefail
cd $TMPDIR
export PATH=/opt/conda/bin:$PATH

/opt/conda/bin/aws s3 cp s3://panorama-clk-repro/${project}/${sample}.Aligned.sortedByCoord.out.bam .
/opt/conda/bin/aws s3 cp s3://panorama-clk-repro/${project}/${sample}.Aligned.sortedByCoord.out.bam.bai .
mkdir GRCh38_star
/opt/conda/bin/aws s3 cp s3://panorama-refs/GRCh38_star/genes.gtf GRCh38_star/


/opt/conda/bin/IsoModule GRCh38_star/genes.gtf ${sample}.Aligned.sortedByCoord.out.bam -i 3 -T 10 -c 1 -v 3 -e 50 -C 20 -o .

/opt/conda/bin/aws s3 cp ${sample}.Aligned.sortedByCoord.out.bam.IsoMatrix s3://panorama-clk-repro/${project}/
/opt/conda/bin/aws s3 cp ${sample}.Aligned.sortedByCoord.out.bam.IsoExon s3://panorama-clk-repro/${project}/

