sudo $(aws ecr get-login --no-include-email --region us-east-1)
sudo docker build -t minimap2 .
sudo docker tag minimap2:latest 977618787841.dkr.ecr.us-east-1.amazonaws.com/minimap2:latest
sudo docker push 977618787841.dkr.ecr.us-east-1.amazonaws.com/minimap2:latest