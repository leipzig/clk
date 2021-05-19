#arn:aws:iam::977618787841:role/service-role/AWSBatchServiceRole
#  "bucket=clk-splicing" "project=SRP091981" \
./batchit submit \
            --image 977618787841.dkr.ecr.us-east-1.amazonaws.com/star:latest \
            --role ecsTaskRole \
            --queue rightnow \
            --jobname biggerSRR5009408 \
            --cpus 16 \
            --mem 32000 \
            --envvars "sample=SRR5009408" "project=SRP091981"\
            --ebs /mnt/my-ebs:500:st1:ext4 \
            star.sh