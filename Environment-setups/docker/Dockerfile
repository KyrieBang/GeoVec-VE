FROM ubuntu:20.04

WORKDIR /home

COPY /home/GeoVec-VE/Environment-setups/dependency /home/dependency

RUN apt-get update &&\
    apt-get install -y vim &&\
    apt-get install -y make &&\
    apt-get install -y cmake &&\
    apt-get install -y g++ &&\
    apt-get install -y python3-pip &&\
    apt-get install -y openjdk-8-jdk &&\
    apt-get install -y maven &&\
    apt-get install -y libmapnik-dev &&\
    apt-get install -y mpich &&\
    pip3 install flask -i https://mirrors.aliyun.com/pypi/simple/ &&\
    pip3 install numpy -i https://mirrors.aliyun.com/pypi/simple/ &&\
    pip3 install pandas -i https://mirrors.aliyun.com/pypi/simple/ &&\
    pip3 install matplotlib -i https://mirrors.aliyun.com/pypi/simple/ &&\
    pip3 install plotly -i https://mirrors.aliyun.com/pypi/simple/ &&\
    pip3 install seaborn -i https://mirrors.aliyun.com/pypi/simple/ &&\
    cd /home/dependency &&\
    tar -zxvf zlib-1.2.11.tar.gz &&\
    cd zlib-1.2.11 && ./configure && make && make install && cd .. && rm -r zlib-1.2.11 &&\
    tar -zxvf libpng-1.2.59.tar.gz &&\
    cd libpng-1.2.59 && ./configure && make && make install && cd .. && rm -r libpng-1.2.59 &&\
    tar -zxvf boost-1.72.0.tar.gz &&\
    cd boost_1_72_0 && ./bootstrap.sh --prefix=/usr/local && ./b2 install --with=all && cd .. && rm -r boost_1_72_0 &&\
    tar -zxvf shapelib-1.5.0.tar.gz &&\
    cd shapelib-1.5.0 && ./configure && make && make install && cd .. && rm -r shapelib-1.5.0 &&\
    tar -zxvf hadoop-2.7.2.tar.gz -C /usr/local &&\
    tar -zxvf spark-3.0.0.tgz -C /usr/local && mv /usr/local/spark-3.0.0-bin-hadoop2.7 /usr/local/spark-3.0.0 &&\
    cp /usr/local/spark-3.0.0/conf/spark-env.sh.template /usr/local/spark-3.0.0/conf/spark-env.sh &&\
    echo 'export SPARK_DIST_CLASSPATH=$(/usr/local/hadoop-2.7.2/bin/hadoop classpath)' >> /usr/local/spark-3.0.0/conf/spark-env.sh &&\
    echo 'export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64' >> ~/.bashrc &&\
    echo 'export HADOOP_HOME=/usr/local/hadoop-2.7.2' >> ~/.bashrc &&\
    echo 'export SPARK_HOME=/usr/local/spark-3.0.0' >> ~/.bashrc &&\
    echo 'export PATH=$HADOOP_HOME/bin:$SPARK_HOME/bin:$PATH' >> ~/.bashrc &&\
    bash -c "source ~/.bashrc" &&\
    cd /home && rm -r /home/dependency


EXPOSE 10085
EXPOSE 8888