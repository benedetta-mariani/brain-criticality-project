import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from avalanches import*
import warnings
warnings.filterwarnings(action = 'ignore')

plt.rcParams['xtick.labelsize']= 20
plt.rcParams['ytick.labelsize']= 20
plt.rcParams['axes.labelsize']= 20
plt.rcParams['legend.fontsize']= 20
plt.rcParams['axes.spines.right'] = True
plt.rcParams['axes.spines.top'] = True
 

## Load Data
data = np.load("data.npy") # data.shape = time x nchannels

## Binarize data
m = np.mean(data,axis = 0)
s = np.std(data,axis = 0)
thres = 3
bindata = threshold(data,m,s,thres,"posneg", "option1")

## Detect avalanches
sizes, durations = main(bindata)

## Plot figure
fig = plt.figure (figsize = (25,6))
fit = pwl.Fit(sizes,xmax = max(sizes), parameter_range= {"alpha":[1,4]},discrete = True)
exp = fit.power_law.alpha
print('Tau is: ', exp)

ax1 = fig.add_subplot(1,3,1)
pwl.plot_pdf(sizes, color = 'tab:blue',marker = 'o')
fit.power_law.plot_pdf(sizes, color = 'tab:blue',linestyle = '--')
ax1.set_xlabel(r'Sizes [$\#$ of events]')
ax1.set_ylabel(r'pdf(sizes)')

ax2 = fig.add_subplot(1,3,2)
fit = pwl.Fit(durations,xmax = max(durations),parameter_range= {"alpha":[1,4]},discrete = True)
exp = fit.power_law.alpha
print('Alpha is: ', value)
pwl.plot_pdf(durations,color = 'tab:red', marker = 'o')
fit.power_law.plot_pdf(durations,color = 'tab:red' ,linestyle = '--')
ax2.set_xlabel(r'Durations [$\#$ of timebins]')
ax2.set_ylabel(r'pdf(durations)')

ax3 = fig.add_subplot(1,3,3)
scaling(sizes,durations,1,4,4,ax3,'default','default','default','default','default','default',max(sizes), max(durations))
ax3.set_xlabel(r'Durations [$\#$ of timebins]')
ax3.set_ylabel(r'Sizes [$\#$ of events]')
plt.savefig('Avalanches.png', dpi = 300,bbox_inches ='tight')
fig.tight_layout()
plt.show()