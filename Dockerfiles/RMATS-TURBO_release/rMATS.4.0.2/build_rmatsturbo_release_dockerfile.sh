sudo $(aws ecr get-login --no-include-email --region us-east-1)
sudo docker build -t rmats-turbo .
sudo docker tag rmats-turbo:latest 977618787841.dkr.ecr.us-east-1.amazonaws.com/rmats-turbo:latest
sudo docker push 977618787841.dkr.ecr.us-east-1.amazonaws.com/rmats-turbo:latest

sudo docker build -t rmats-turbo .
sudo docker tag rmats-turbo:latest quay.io/leipzig/rmats-turbo:latest
sudo docker push quay.io/leipzig/rmats-turbo