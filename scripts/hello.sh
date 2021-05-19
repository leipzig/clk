set -euo pipefail
cd $TMPDIR

echo "hello world" > helloworld.txt
aws s3 cp helloworld.txt s3://clk-splicing/helloworld.txt
