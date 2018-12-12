tar -c -f v2.0.3.tar rmats2sashimiplot-2.0.3

sudo $(aws ecr get-login --no-include-email --region us-east-1)
sudo docker build -t rmats-sashimi .
sudo docker tag rmats-sashimi:latest 977618787841.dkr.ecr.us-east-1.amazonaws.com/rmats-sashimi:latest
sudo docker push 977618787841.dkr.ecr.us-east-1.amazonaws.com/rmats-sashimi:latest

sudo docker build -t rmats-sashimi .
sudo docker tag rmats-sashimi:latest quay.io/leipzig/rmats-sashimi:latest
sudo docker push quay.io/leipzig/rmats-sashimi