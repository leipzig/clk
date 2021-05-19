./batchit submit \
            --image worker:latest \
            --role AWSBatchServiceRole \
            --queue queue_leipzig \
            --jobname SRR5009427 \
            --cpus 32 \
            --mem 32000 \
            --envvars "sample=SRR5009427" "bucket=clk-splicing/SRP091981" \
            --ebs /mnt/my-ebs:500:st1:ext4 \
            star.sh