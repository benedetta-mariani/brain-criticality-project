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
is a C++ (ROOT) implementation of ```nunits``` decoupled Ornstein Uhlenbeck processes with a common diffusion coefficient that varies in time (here accordingly to a thresholded Ornstein Uhlenbeck process), as a minimal model for neural activity of independent neurons with a common external driver. Each process evovles according to the equation:
$ x_i(t + \Delta t) = x_i(t) - \frac{x_i}{tau}dt + \sqrt{D(t)}dW(t)$


After compliling, it is obviously way much faster than the Python implementation.

To compile it and run it:

```
root -l
.L OrnsteinUhlenbeck.C++
int ntime = 100;
int nunits = 200;
Main(nunits,ntime)
```
To open the results in Python and analyze them just run:

```
import numpy as np
nunits = 200
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
            neuro[i] = np.array(row)
            i = i +1
            
```






