sudo $(aws ecr get-login --no-include-email --region us-east-1)
sudo docker build -t star .
sudo docker tag star:latest 977618787841.dkr.ecr.us-east-1.amazonaws.com/star:latest
sudo docker push 977618787841.dkr.ecr.us-east-1.amazonaws.com/star:latest