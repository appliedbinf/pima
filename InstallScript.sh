#!/usr/bin/env bash

# Install all the Docker Environment
sudo apt-get update

sudo apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common

wget -qO- https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

wget https://raw.githubusercontent.com/appliedbinf/pima-docker2/main/pima_interface.py

sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io

# Install the Nvidia Toolkit and Nvidia Drivers
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)

wget -qO- https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
wget -qO- https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | sudo tee /etc/apt/sources.list.d/nvidia-docker.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit

docker pull appliedbioinformaticslab/pimadocker2:latest
