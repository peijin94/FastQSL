# FastQFact

A module to calculate Q Factor with GPU.

The idea is to do the most computational intensive work in GPU with compiled code (TraceBline.cu), and do the rest of the complex but not computational intensive work in Python.

## Dependencies

* CUDA >= 10.2
* PyCuda
* Cupy

## Step by step Tutorial

### Preparation

#### Hardware

Buy a computer with Nvidia GPU

#### Software

* Install GPU driver, 驱动精灵[http://www.drivergenius.com/](http://www.drivergenius.com/) is a easy way for windows, for 'apt' based Linux, try ```sudo ubuntu-drivers autoinstall```.
* [https://developer.nvidia.com/Cuda-downloads](https://developer.nvidia.com/Cuda-downloads) Download cuda **v10.2! (don't install the newest)**
* Install anaconda [https://docs.anaconda.com/anaconda/install/windows/](https://docs.anaconda.com/anaconda/install/windows/)
* Install git
* Install c/c++ compiler, for Windows [https://visualstudio.microsoft.com/visual-cpp-build-tools](https://visualstudio.microsoft.com/visual-cpp-build-tools), for linux ```sudo apt install gcc```

#### Build Enviroment

##### Test conda

Start a 'Anaconda Powershell prompt' (in Windows), a terminal (in Linux). Type command:

```bash
conda --version
```

You should get "conda x.x.x", which means conda is successfully installed.

##### New conda-env

Create a new env named 'fastqsl' and use 3.8.2 version of python:

```bash
conda create -n fastqsl python=3.8.2
conda activate fastqsl
```

after this you should be able to see a '(fastqsl)' in your command line.

##### Install dependencies

```bash
conda install numpy matplotlib pyvista scipy jupyterlab tqdm
python -m pip install cupy-cuda102 ipython-autotime 
```

### Run

#### Download code

```bash
git clone https://github.com/Pjer-zhang/FastQSL.git
```

#### Run in jupyterlab

In the env created last section, type command

```bash
cd FastQSL
python -m jupyterlab --port 9999
```

Normally it will pop up a browser automatically, if it didn't, you need to it manually by copy the url in the output of the command to a browser (Chrome recomended)

Enjoy.
