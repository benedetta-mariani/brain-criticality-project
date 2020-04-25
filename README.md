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

- to detect avalanches and plot their sizes and lifetimes distributions simply run

```
sizes, lifetimes = main(dataset,interv)
plt.figure()
pwl.plot_pdf(sizes, color='r', linewidth=2, label='pdf', linear_bins = False)
plt.figure()
pwl.plot_pdf(lifetimes, color='b', linewidth=2, label='pdf', linear_bins = False)
```
   where ```dataset``` in a tri/bidimensional arrya (temporal dim x spatial dim) of discretized signal and ```interv``` is the width of the chosen temporal bin to bin the data.
   
- Detailed instructions in the tutorial notebook
