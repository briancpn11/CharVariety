# CharVariety

## Introduction

CharVariety is a complementary library to the following paper on arxiv:

[Stationary measures and orbit closures of uniformly expanding random dynamical systems on surfaces](https://arxiv.org/pdf/2006.03166.pdf) by <cite> Ping Ngai (Brian) Chung </cite>

The goal is to implement the algorithm in Chapter 8 which allows one to verify <u>uniform expansion</u> on grid points of a given character variety (see Chapter 8 of the paper for details). 

See Theorem C and D for key theoretical consequences of <u>uniform expansion</u> and motivation. 

</br></br></br>

## Versions

Currently there are three versions of the code:

1. <b>Basic version</b>: `main/mainStandAloneScript.cpp`

This can be compiled without any extra package other than basic g++. To compile on terminal: at the CharVariety directory, use

```bash
g++ -o .\main\mainStandAloneScript .\main\mainStandAloneScript.cpp
```

2. <b>With multithreading</b>: `main/mainStandAloneWithMultithread.cpp`

This is the same as the basic version except it allows for multithreading using OpenMP. To compile:

```bash
g++ -o .\main\mainStandAloneWithMultithread .\main\mainStandAloneWithMultithread.cpp -fopenmp
```

3. <b>Latest Version (as of September 13th, 2021) </b>: `main/mainScript.cpp`

I am in the process of refactoring the script to make it more readable and potentially faster. This is currently the latest version. To compile:

```bash
g++ -g .\main\mainScript.cpp .\utils\helpers.cpp -o .\main\mainScript -fopenmp
```

4. <b>Interval Arithmetic</b>: to come

To make sure the result is not affected by floating point error, all the arithmetic should be done using interval arithmetic. Currently the plan is to use the `boost/numeric/interval` package. 

</br></br></br>

## Results

No matter which version is compiled, the output is stored in the `result` directory as a txt file. It should look like the following:

```bash
Program run on: Mon Sep 13 16:11:22 2021
Minimum average expansion: 0.531619
Other parameters: 
 k-2: 1.990000
 Spacing in the surface directions: r: 0.100000
 Spacing in the tangent direction: rr: 0.010000
 Total number of grid points on the surface: 5190
The point where the min average expansion occur: P = (1.928409, 1.900000, 1.700000)
Time taken: 1.879000 seconds
```
The ''Minimum average expansion'' is positive means that uniform expansion is verified on the grid points. 

## Unit tests
There are unit tests in the `unittest` directory, implemented using Boost.Test (I am using Boost Version 1.75.0).