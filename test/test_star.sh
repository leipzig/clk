aws s3 cp s3://clk-splicing/SRP091981/${sample}_1.fastq.gz . 
aws s3 cp s3://clk-splicing/SRP091981/${sample}_2.fastq.gz .
STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --alignEndsType EndToEnd --genomeDir GRCh38_star --runThreadN ${cpus} --outFileNamePrefix ${sample}.  --readFilesIn  ${sample}_1.fastq.gz ${sample}_2.fastq.gz
