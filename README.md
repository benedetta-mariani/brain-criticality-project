# Avalanches.py 
is a short module to detect events in time series of neural activity, plot raster plots and  calculate avalanches statistics.


## Setup:

- get the module

Packages needed:

- numpy

- matplotlib.pyplot

- powerlaw (https://github.com/jeffalstott/powerlaw)


Other functions needed:
- scipy.stats.percentileofscore

- scipy.optimize.curve_fit

- scipy.signal.find_peaks

- statsmodels.regression.linear_model

## to use it: 


```from avalanches import * ```

- to detect avalanches and calculate their sizes and lifetimes simply run:

```
sizes, lifetimes = main(dataset,interv)
```
   where ```dataset``` is a tri/bidimensional array (shape = temporal dim x spatial dims) of discretized signals and ```interv``` is the width of the chosen temporal bin to bin the data.
   
- Detailed instructions in the tutorial notebook

# OrnsteinUhlenbeck.C
is a C++ (ROOT) implementation of N Ornstein Uhlebenck processes with a diffusion coefficient that varies in time (here accordingly to another Ornstein Uhlenbeck porcess thresholded). After compliling, it is obviously much faster than the Python implementation.

To compile it and run it:

```
root -l
.L OrnsteinUhlenbeck.C++
int ntime = 1000;
int nunits = 200;
Main(nunits,ntime)
```
To open the results in Python and analyze them just run:

```
nunits = 200
ntime = 1000
dt = 0.001
N = int(ntime/dt)
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






