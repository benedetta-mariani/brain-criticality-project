# Avalanches.py 
is a short module to detect events in time series of neural activity, plot raster plots and  calculate avalanches statistics.


## Setup:

- get the repository:

```git clone https://github.com/benedetta-mariani/Avalanches-in-rat-barrel-cortex```
- to use it: 
```from avalanches import * ```


It needs:

- numpy

- matplotlib.pyplot

- powerlaw (https://github.com/jeffalstott/powerlaw)


Functions needed:
- scipy.stats.percentileofscore

- scipy.optimize.curve_fit

- scipy.signal.find_peaks
