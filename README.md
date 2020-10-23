# Discrete Optimal Transport Library

[![Build Status](https://travis-ci.org/stegua/dotlib.svg?branch=master)](https://travis-ci.org/stegua/dotlib)

This repository contains basic scripts for computing Kantorovich.Wasserstein distances.
The long term objective is to develop Discrete Optimal Transport Algorithms.

# Notebooks

| Data | Notebook | Link |
|:-|:-|:-|
|**[2020/10/23]**|*Two LP models for computing KW-distances*|[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dotlib/.ipynb)|

### Early stage of the project
The project started in 2016 and it containes old python scripts to compute the distances 
between pair of images of the [DOTMark benchmarks](http://www.stochastik.math.uni-goettingen.de/index.php?id=215/).

To run the Python scripts you need the following libraries:

* [Matplotlib](https://matplotlib.org/) (python)
* [NetworkX](http://networkx.github.io/) (python)
* [Gurobi](http://www.gurobi.com)

To build the library from C++ source, you need the following external libraries:

* [Lemon Graph Library](http://lemon.cs.elte.hu/trac/lemon)

## Reference
[1] Bassetti F., Gualandi S., Veneroni M. *On the Computation of Kantorovich-Wasserstein Distances between 2D-Histograms by Uncapacitated Minimum Cost Flows*. [Available on arXiv](https://arxiv.org/abs/1804.00445). Submitted on April, 2nd,  2018.


### MIT License
Copyright (c) 2017-2020, by Stefano Gualandi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
