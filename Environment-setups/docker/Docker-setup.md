
## Installation instructions for running Docker image on Ubuntu


### 1. Install Docker
* [Docker](https://www.docker.com/) ( recommended version 27.2.0 )
> ~~~
> sudo apt-get update
> sudo apt install apt-transport-https ca-certificates curl software-properties-common gnupg lsb-release
> curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
> sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
> sudo apt-get install docker-ce docker-ce-cli containerd.i
> sudo docker version
> ~~~

### 2. Download the Docker image
> ~~~
> docker pull crpi-dvwy5e9a3vrsg1gg.cn-guangzhou.personal.cr.aliyuncs.com/kyrie-bang/geove-ve-env:v1.0
> docker images
> ~~~



### 3. Create the Docker container
> ~~~
> mkdir /home/GeoVec-VE
> docker run -it --name geovec-ve-env -p 8888:8888 -p 10085:10085 -v /home/GeoVec-VE/:/home/GeoVec-VE geove-ve-env:v1.0
> (Docker) cd ./GeoVec-VE
> ~~~