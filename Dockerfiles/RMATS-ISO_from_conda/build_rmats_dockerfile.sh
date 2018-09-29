sudo $(aws ecr get-login --no-include-email --region us-east-1)
sudo docker build -t rmats-iso .
sudo docker tag rmats-iso:latest 977618787841.dkr.ecr.us-east-1.amazonaws.com/rmats-iso:latest
sudo docker push 977618787841.dkr.ecr.us-east-1.amazonaws.com/rmats-iso:latest