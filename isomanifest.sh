set -euo pipefail
cd $TMPDIR

echo "downloading bams.."
for f in ${untreated}; do aws s3 cp s3://$f .; done;
for f in ${treated}; do aws s3 cp s3://$f .; done;

aws s3 cp s3://${manifest} manifest.list

echo "I am here:"
pwd
ls -alt .

echo "downloading gtf.."
mkdir GRCh38_star
aws s3 cp s3://panorama-refs/GRCh38_star/${gtf} GRCh38_star/

python3 /rMATS-ISO-master/rMATS-ISO.py  --in-gtf GRCh38_star/${gtf} --in-bam manifest.list -o ${comparison}

aws s3 sync ${comparison} s3://panorama-clk-repro/${project}/${comparison}

