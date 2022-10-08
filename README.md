# FastQSL

A module to calculate Q Factor with GPU, [FastQSL](https://arxiv.org/abs/2208.12569)

The idea is to do the most computational intensive work in GPU with compiled code (TraceBline.cu), and do the rest of the complex but not computational intensive work in Python.

A CPU version is also provided [https://github.com/el2718/FastQSL](https://github.com/el2718/FastQSL).

## Dependencies

* CUDA >= 11
* PyCuda
* Cupy

## Step by step Tutorial

### Preparation

#### Hardware

Buy a computer with Nvidia GPU

#### Software

* Install GPU driver, [https://www.nvidia.com/Download/index.aspx](https://www.nvidia.com/Download/index.aspx), for 'apt' based Linux, try ```sudo ubuntu-drivers autoinstall```.
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


-----------------------------

## Cite as

* Zhang, P., Chen, J.*, Liu, R. and Wang, C., 2022, FastQSL: A Fast Computation Method for Quasi-separatrix Layers. The Astrophysical Journal, 937, 26

```bibtex
@ARTICLE{2022ApJ...937...26Z,
       author = {{Zhang}, PeiJin and {Chen}, Jun and {Liu}, Rui and {Wang}, ChuanBing},
        title = "{FastQSL: A Fast Computation Method for Quasi-separatrix Layers}",
      journal = {\apj},
     keywords = {Solar magnetic fields, GPU computing, 1503, 1969},
         year = 2022,
        month = sep,
       volume = {937},
       number = {1},
          eid = {26},
        pages = {26},
          doi = {10.3847/1538-4357/ac8d61},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022ApJ...937...26Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
