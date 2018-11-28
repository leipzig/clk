set -euo pipefail
cd $TMPDIR

echo "hello world" > helloworld.txt
aws s3 cp helloworld.txt s3://panorama-clk-repro/helloworld.txt
