# Avalanches.py 
is a short module to detect events in time series of neural activity, plot raster plots and  calculate avalanches statistics.


## Setup:

- get the module

Main packages needed:

- numpy

- matplotlib.pyplot

- powerlaw (https://github.com/jeffalstott/powerlaw)

- some functions of scipy

## to use it: 


```from avalanches import * ```

- to detect avalanches and calculate their sizes and lifetimes simply run:

```
sizes, lifetimes = main(dataset,interv)
```
   where ```dataset``` is a tri/bidimensional array (shape = temporal dim x spatial dims) of discretized signals and ```interv``` is the width of the chosen temporal bin to bin the data.
   
- Detailed instructions in the tutorial notebook

# OrnsteinUhlenbeck.C
is a C++ (ROOT) implementation of a multivariate Ornstein Uhlenbeck process on a network with a diffusion coefficient that varies in time (here accordingly to a thresholded Ornstein Uhlenbeck process), as a minimal model for neural activity. The process at node i evolves according to the equation:

![equation](http://www.sciweavers.org/download/Tex2Img_1601221354.jpg)

where a<sub>ij</sub> is the synaptic strength between unit i and j, i.e. the element of the adjacency matrix. 

Note that this equation is entirely equivalent to the linearized version of a noisy neural network of Wilson-Cowan type [[1]](https://seis.bristol.ac.uk/~sb15704/papers/267384.pdf).

In the case in which every a<sub>ij</sub> is set zero the model reduces to ```nunits``` decoupled Ornstein Uhlenbeck processes with a commmon time varying diffusion coefficient.

Note that the spectral radius of the adjacency matrix must be less than 1, in order not to have signals that grow/decay in time without bounds [[2]](https://arxiv.org/pdf/2007.07447.pdf). 

After compliling, it is obviously way much faster than the Python implementation.

To compile it and run it:

```
root -l
.L OrnsteinUhlenbeck.C++
int ntime = 100;
int nunits = 220;
Main(nunits,ntime)
```
To open the results in Python and analyze them just run:

```
import numpy as np
nunits = 220
ntime = 100
dt = 0.001
tau = 0.08
N = int(ntime/dt) - int(100*tau/dt)
neuro = np.array([[0. for r in range(N)] for s in range(nunits)])

import csv
with open('output.csv') as csv_file:
    
    csv_reader = csv.reader(csv_file, delimiter=',')

    i = 0
    for row in csv_reader:
        if i < nunits:
            neuro[i] = np.array(row[int(100*tau/dt):]) #to avoid initial transient
            i = i +1
            
```






