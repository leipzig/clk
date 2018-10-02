set -euo pipefail
cd $TMPDIR

#for some reason aws batch doesn't get this path
export PATH=$PATH:/opt/conda/bin

mkdir GRCh38_star
echo "download genes..."
aws s3 cp s3://panorama-refs/GRCh38_star/genes.gtf GRCh38_star/genes.gtf
echo "download sam..."
aws s3 cp s3://panorama-clk-repro/${project}/${sample}.sam .
echo "download rRNA gtf"
aws s3 cp s3://panorama-refs/GRCh38_star/rRNA_tx.gtf GRCh38_star/rRNA_tx.gtf 

# "{params.lr2rmats} filter {input.sam} {params.rm_gtf} -v {params.aln_cov} -q {params.iden_frac} -s {params.sec_rat} 2> {log} | {params.samtools} sort -@ {threads} > {output.filtered_bam} 2>> {log}; "
# "{params.lr2rmats} update-gtf {output.filtered_bam} {input.gtf} 2>> {log} > {output.sam_gtf}"
echo "lr2rmats filter..."    
#lr2rmats filter {sample}.sam GRCh38_star/rRNA.gtf -v {aln_cov} -q {iden_frac} -s {sec_rat} 2> {log} | samtools sort -@ {cpus} > {sample}.filtered_bam 2>> {sample}.lr2rmats.log
echo "lr2rmats filter {sample}.sam GRCh38_star/rRNA.gtf -v {aln_cov} -q {iden_frac} -s {sec_rat}  | samtools sort -@ {cpus} > {sample}.filtered_bam"
#lr2rmats filter ${sample}.sam -r GRCh38_star/rRNA.gtf -v ${aln_cov} -q ${iden_frac} -s ${sec_rat}  | samtools sort -@ ${cpus} > {sample}.filtered_bam

#[trans_exon_comp] Strands of exons do NOT match.
#lr2rmats filter ${sample}.sam -r GRCh38_star/rRNA.gtf -v ${aln_cov} -q ${iden_frac} -s ${sec_rat}  | samtools sort -@ ${cpus} > {sample}.filtered_bam

#lets use a rRNA.gtf with no negative strnads
lr2rmats filter ${sample}.sam -r GRCh38_star/rRNA_tx.gtf -v ${aln_cov} -q ${iden_frac} -s ${sec_rat}  | samtools sort -@ ${cpus} > ${sample}.filtered_bam

echo "yo" > ${sample}.lr2rmats.log
echo "lr2rmats update-gtf..."  
lr2rmats update-gtf ${sample}.filtered_bam GRCh38_star/genes.gtf > ${sample}_sam_novel.gtf
#2>> ${sample}.lr2rmats.log
aws s3 cp ${sample}_sam_novel.gtf s3://panorama-clk-repro/${project}/
aws s3 cp ${sample}.filtered_bam s3://panorama-clk-repro/${project}/
aws s3 cp ${sample}.lr2rmats.log s3://panorama-clk-repro/${project}/