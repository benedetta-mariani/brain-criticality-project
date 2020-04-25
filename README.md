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
fig = plt.figure(figsize = (13,4))
ax1 = fig.add_subplot(1,2,1)
pwl.plot_pdf(sizes, color='r', linewidth=2, label='pdf', linear_bins = False,ax = ax1)
ax2 = fig.add_subplot(1,2,2)
pwl.plot_pdf(durations, color='b', linewidth=2, label='pdf', linear_bins = False,ax = ax2)
```
   where ```dataset``` is a tri/bidimensional array (shape = temporal dim x spatial dims) of discretized signals and ```interv``` is the width of the chosen temporal bin to bin the data.
   
- Detailed instructions in the tutorial notebook
