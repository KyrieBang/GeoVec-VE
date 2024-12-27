
## Installation Instructions to deploy environment on Ubuntu

### 1. base environments
> ~~~
> sudo apt install vim
> sudo apt install make
> sudo apt install cmake
> sudo apt install g++
> sudo apt-get install python3-pip 
> sudo apt install openjdk-8-jdk
> sudo apt install maven
> echo ‘export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64’ >> ~/.bashrc
> source ~/.bashrc
> ~~~

### 2. C++
* [MPICH](https://www.mpich.org/) ( recommended version 3.3.2 )
> ~~~
> wget http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz
> tar -zxvf mpich-3.3.2.tar.gz
> cd mpich-3.3.2
> ./configure --disable-fortran 
> make && make install
> ~~~

* [zlib](http://www.zlib.net/) ( recommended version 1.2.11 )
> ~~~
> wget http://www.zlib.net/zlib-1.2.11.tar.gz
> tar -zxvf zlib-1.2.11.tar.gz
> cd zlib-1.2.11
> ./configure
> make && make install
> ~~~
  
* [libpng](http://www.libpng.org/pub/png//libpng.html) ( recommended version 1.2.59 )
> ~~~
> wget https://sourceforge.net/projects/libpng/files/libpng12/1.2.59/libpng-1.2.59.tar.gz
> tar -zxvf libpng-1.2.59.tar.gz
> cd libpng-1.2.59
> ./configure
> make && make install
> ~~~

* [mapnik](https://mapnik.org) ( recommended version 3.0.23 )
> ~~~
> sudo apt-get install libmapnik-dev
> ~~~

* [Boost](https://www.boost.org) ( recommended version 1.7.2 )
> ~~~
> wget https://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.gz
> tar -zxvf boost_1_72_0.tar.gz
> cd boost_1_72_0
> ./bootstrap.sh --prefix=/usr/local
> ./b2 install --with=all
> ~~~

* [shapelib](http://shapelib.maptools.org/) ( recommended version 1.5.0 )
> ~~~
> wget https://download.osgeo.org/shapelib/shapelib-1.5.0.tar.gz
> tar -zxvf shapelib-1.5.0.tar.gz
> cd shapelib-1.5.0
> ./configure
> make && make install
> ~~~

* [Crow](https://github.com/ipkn/crow) ( recommended version 0.1 )

    header files are saved in "./Algorithms/include/dep-include/crow".

* [Hicore](https://gitee.com/CoreSpatial/core-open) ( closed-source )

    header files and binary files are saved in "./Algorithms/include/dep-include/hicore" and "./Algorithms/library/dep-lib".


### 3. Python

* [flask](https://flask.palletsprojects.com/en/2.0.x/) ( recommended version 2.0.3 )
> ~~~
> pip install flask -i https://mirrors.aliyun.com/pypi/simple/
> ~~~

* [numpy](https://numpy.org/) ( recommended version 1.24.4 )
> ~~~
> pip install numpy -i https://mirrors.aliyun.com/pypi/simple/
> ~~~

* [pandas](https://pandas.pydata.org/) ( recommended version 2.0.3 )
> ~~~
> pip install pandas -i https://mirrors.aliyun.com/pypi/simple/
> ~~~

* [matplotlib](https://matplotlib.org/) ( recommended version 3.7.5 )
> ~~~
> pip install matplotlib -i https://mirrors.aliyun.com/pypi/simple/
> ~~~

* [plotly](https://plotly.com/python/) ( recommended version 5.23.0 )
> ~~~
> pip install plotly -i https://mirrors.aliyun.com/pypi/simple/
> ~~~

* [seaborn](https://seaborn.pydata.org/) ( recommended version 0.13.2 )
> ~~~
> pip install seaborn -i https://mirrors.aliyun.com/pypi/simple/
> ~~~


### 4. JavaScripts

* [openLayers](https://openlayers.org/) ( recommended version 6.13.0 )
> ~~~
> npm install ol
> ~~~

* [Echarts](https://echarts.apache.org/zh/index.html) ( recommended version 5.5.1 )
> ~~~
> npm install echarts
> ~~~


### 5. Java

* [Hadoop](https://hadoop.apache.org/) ( recommended version 2.7.2 )
> ~~~
> wget https://archive.apache.org/dist/hadoop/core/hadoop-2.7.2/hadoop-2.7.2.tar.gz
> tar -zxvf hadoop-2.7.2.tar.gz -C /usr/local
> echo ‘export HADOOP_HOME=/usr/local/hadoop-2.7.2’ >> ~/.bashrc
> echo ‘export PATH=$HADOOP_HOME/bin:$PATH’ >> ~/.bashrc
> source ~/.bashrc
> ~~~

* [Spark](https://spark.apache.org/) ( recommended version 3.0.0 )
> ~~~
> wget https://archive.apache.org/dist/spark/spark-3.0.0/spark-3.0.0-bin-hadoop2.7.tgz
> tar -zxvf spark-3.0.0-bin-hadoop2.7.tgz -C /usr/local
> mv /usr/local/spark-3.0.0-bin-hadoop2.7 /usr/local/spark-3.0.0
> cp /usr/local/spark-3.0.0/conf/spark-env.sh.template /usr/local/spark-3.0.0/conf/spark-env.sh
> echo ‘export SPARK_DIST_CLASSPATH=$(/usr/local/hadoop-2.7.2/bin/hadoop classpath)’ >> /usr/local/spark-3.0.0/conf/spark-env.sh
> echo ‘export SPARK_HOME=/usr/local/spark-3.0.0’ >> ~/.bashrc
> echo ‘export PATH=$SPARK_HOME/bin:$PATH’ >> ~/.bashrc
> source ~/.bashrc
> ~~~