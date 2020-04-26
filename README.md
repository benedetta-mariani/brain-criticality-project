# Avalanches.py 
is a short module to detect events in time series of neural activity, plot raster plots and  calculate avalanches statistics.


## Setup:

- get the module

It needs:

- numpy

- matplotlib.pyplot

- powerlaw (https://github.com/jeffalstott/powerlaw)


Functions needed:
- scipy.stats.percentileofscore

- scipy.optimize.curve_fit

- scipy.signal.find_peaks

## to use it: 


```from avalanches import * ```

- to detect avalanches and calculate their sizes and lifetimes simply run:

```
sizes, lifetimes = main(dataset,interv)
```
   where ```dataset``` is a tri/bidimensional array (shape = temporal dim x spatial dims) of discretized signals and ```interv``` is the width of the chosen temporal bin to bin the data.
   
- Detailed instructions in the tutorial notebook
