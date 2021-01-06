import numpy as np
import matplotlib.pyplot as plt

from avalanches import*


##Data:

lista = [r"/home/bm/Desktop/altriDati/20181018_4x64/2s/ISI_2s_1st_sf976.5625_001.bin"]

## Load Data

path = lista[0]
width = 4    #columns
heigth = 64   #rows
n_words_to_load = -1  #load  all frames

raw = np.fromfile( path, dtype=np.double, count = n_words_to_load)  
a = np.reshape(raw, (-1, width, heigth))

a1 = []
for z in range(len(a)):
	a1.append(np.delete(a[z],(0,1,2,3,4,5,6,7,8), axis = 1).tolist())

a1 = np.asarray(a1)
nunits = width*(heigth-9)
a1 = a1.reshape(-1,nunits)

## Binarize them

bindata = threshold(a1,np.mean(a1,axis = 0),np.std(a1,ddof = 1,axis = 0),3,"posneg", "option1")

## Detect avalanches

sizes, durations = main(bindata)

draw(sizes, max(sizes), max(sizes), "power_law")
plt.xlabel ('Sizes [# of events]')
plt.ylabel('pdf(sizes)')
plt.show()

draw(durations, max(durations), max(durations),"power_law")
plt.xlabel ('Durations [# of time bins]')
plt.ylabel('pdf(durations)')
plt.show()

scaling(sizes,durations,1)
plt.show()

