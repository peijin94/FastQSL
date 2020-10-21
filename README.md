# FastQFact

A module to calculate Q Factor with GPU.

The idea is to do the most computational intensive work in GPU with compiled code (TraceBline.cu), and do the rest of the complex but not computational intensive work in Python.

## Dependency

* CUDA >= 10.2
* PyCuda
