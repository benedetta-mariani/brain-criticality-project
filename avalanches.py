import numpy as np
import matplotlib.pyplot as plt
import powerlaw as pwl
from scipy.stats import percentileofscore
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


def threshold(sample1,means,stds,thres,choose):
    """ 
    Detects as events the points of maximum excursion over a threshold, considering either positive and negative excursions or only negative. Differs from "Fasterthreshold" since here only the one largest maximum between two crossings of the mean assigns the final event time.
    
    Parameters
    --------
    sample1 : tridimensional array (shape = temporal dim x spatial dim1 x  spatial dim2 ) of recorded voltages
    means : bidimensional array (shape = spatial dim1 x spatial dim2 ) of the means of the signals
    stds : bidimensional array (shape = spatial dim1 x spatial dim2 ) of the thresholds for each channel
    thres : multiplicative coefficent for the thresholds
    choose : if "posneg" both positive and negative deflections are considered, if "neg" only negative
    
    Returns
    --------
    sample2 : discretized tridimensional array 
    """
    
    if sample1.ndim < 3:
        sample1 = sample1.reshape(sample1.shape[0],sample1.shape[1],1)
        means = means.reshape(sample1.shape[1],1)
        stds = stds.reshape(sample1.shape[1],1)
    
    sample2 = np.zeros(sample1.shape, dtype = int)
    

    if choose == "posneg":
        

        for i in range(sample1.shape[1]):
            for j in range(sample1.shape[2]):
                if stds[i][j] > 0:
                    picchi = [[]]
                    tempi = [[]]
                    n = 0
                    for z in range(sample1.shape[0]):
                            if np.abs(sample1[z,i,j]- means[i][j]) >=  thres*stds[i][j]:
                                picchi[n].append(np.abs(sample1[z,i,j]- means[i][j]))
                                tempi[n].append(z)

                            else:
                                if z >0:
                                    if np.sign(sample1[z,i,j]- means[i][j]) !=  np.sign(sample1[z-1,i,j]- means[i][j]): 

                                        n = n + 1
                                        picchi.append([])
                                        tempi.append([])

                    zeta = []       
                    for k in range(picchi.count([])):
                        picchi.remove([])
                    for p in range(tempi.count([])):
                        tempi.remove([])
                    for l in range(len(picchi)):
                        zeta.append(tempi[l][picchi[l].index(max(picchi[l]))])
                    for r in range(len(sample1)):
                        if r in zeta:
                            sample2[r,i,j] = 1
                        else:
                            sample2[r,i,j] = 0         
                        
        
    if choose == "neg":
    
        for i in range(sample1.shape[1]):
            for j in range(sample1.shape[2]):
                if stds[i][j] > 0:
                    picchi = [[]]
                    tempi = [[]]
                    n = 0
                    for z in range(sample1.shape[0]):
                        if (sample1[z,i,j]- means[i][j]) <= - thres*stds[i][j]:
                                picchi[n].append(np.abs(sample1[z,i,j]- means[i][j]))
                                tempi[n].append(z)

                      
                        else:
                            if z > 0:
                                if np.sign(sample1[z ,i,j]- means[i][j]) !=  np.sign(sample1[z-1,i,j]- means[i][j]):
                                    n = n + 1
                                    picchi.append([])
                                    tempi.append([])

                    zeta = []       
                    for k in range(picchi.count([])):
                        picchi.remove([])
                    for p in range(tempi.count([])):
                        tempi.remove([])
                    for l in range(len(picchi)):
                     
                        zeta.append(tempi[l][picchi[l].index(max(picchi[l]))])

                    for r in range(sample1.shape[0]):
                        if r in zeta:
                            sample2[r,i,j] = 1
                        else:
                            sample2[r,i,j] = 0         

            
    return sample2

def fasterthreshold(sample1,means,stds,thres,ref,choose = "posneg"):
    """ 
    Detects as events the points of maximum excursion over a threshold, considering either positive and negative excursions or only negative. Works without exceptions, but still is slow. For a really fast thresholding use findpeaks
    
    Parameters
    --------
    sample1 : tridimensional array (shape = temporal dim x spatial dim1 x  spatial dim2 ) of recorded voltages
    means : bidimensional array (shape = spatial dim1 x spatial dim2 ) of the means of the signals
    stds : bidimensional array (shape = spatial dim1 x spatial dim2 ) of the thresholds for each channel
    thres : multiplicative coefficent for the thresholds
    ref : refractory period. Minimum distance between consecutive events.
    choose : if "posneg" both positive and negative deflections are considered, if "neg" only negative
    
    Returns
    --------
    sample2 : discretized tridimensional array 
    """
    if sample1.ndim < 3:
        sample1 = sample1.reshape(sample1.shape[0],sample1.shape[1],1)
        means = means.reshape(sample1.shape[1],1)
        stds = stds.reshape(sample1.shape[1],1)
        
    sample2 = np.zeros(sample1.shape, dtype = int)
    
    if choose == "posneg":

        for i in range(sample1.shape[1]):
            for j in range(sample1.shape[2]):
                picchi = [[]]
                tempi = [[]]
                n = 0
                for z in range(sample1.shape[0]):
                        if np.abs(sample1[z,i,j]- means[i][j]) >=  thres*stds[i][j] and stds[i][j] > 0:
                            picchi[n].append(np.abs(sample1[z,i,j]- means[i][j]))
                            tempi[n].append(z)
                            if z + 1 < len(sample1):

                                if np.abs(sample1[z + 1,i,j]- means[i][j]) < thres*stds[i][j]:
                                    n = n + 1
                                    picchi.append([])
                                    tempi.append([])

                zeta = []       
                for k in range(picchi.count([])):
                    picchi.remove([])
                for p in range(tempi.count([])):
                    tempi.remove([])
                for l in range(len(picchi)):
                    if l != 0:
                        if tempi[l][picchi[l].index(max(picchi[l]))] - tempi[l-1][picchi[l-1].index(max(picchi[l-1]))] > ref:
                            zeta.append(tempi[l][picchi[l].index(max(picchi[l]))])
                    else:
                        zeta.append(tempi[l][picchi[l].index(max(picchi[l]))])
                for r in range(sample1.shape[0]):
                    if r in zeta:
                        sample2[r,i,j] = 1
                    else:
                        sample2[r,i,j] = 0  
                        
    if choose == "neg":
    
        for i in range(sample1.shape[1]):
            for j in range(sample1.shape[2]):
                picchi = [[]]
                tempi = [[]]
                n = 0
                for z in range(sample1.shape[0]):
                        if (sample1[z,i,j]- means[i][j]) <= - thres*stds[i][j] and stds[i][j] > 0:
                            picchi[n].append(np.abs(sample1[z,i,j]))
                            tempi[n].append(z)
                            if z + 1 < len(sample1):

                                if (sample1[z + 1,i,j]- means[i][j]) > - thres*stds[i][j]: #check
                                    n = n + 1
                                    picchi.append([])
                                    tempi.append([])

                zeta = []       
                for k in range(picchi.count([])):
                    picchi.remove([])
                for p in range(tempi.count([])):
                    tempi.remove([])
                for l in range(len(picchi)):
                    if l != 0:
                        if tempi[l][picchi[l].index(max(picchi[l]))] - tempi[l-1][picchi[l-1].index(max(picchi[l-1]))] > ref:
                            zeta.append(tempi[l][picchi[l].index(max(picchi[l]))])
                    else:
                        zeta.append(tempi[l][picchi[l].index(max(picchi[l]))])
                for r in range(sample1.shape[0]):
                    if r in zeta:
                        sample2[r,i,j] = 1
                    else:
                        sample2[r,i,j] = 0         
    
    return sample2

def findpeaks(sig, thres, choose):
    """
    Finds peaks in a time series.
    Set eventual other constraints in the scipy.signal function find_peaks
    """
    sig2 = np.zeros(sig.shape, dtype = int)
    if choose == "neg":
        p1 = find_peaks(-(sig - np.mean(sig)), height = thres)[0]
        if np.any(p1):
            sig2[p1] = 1
            
    if choose == "posneg":
        p1 = find_peaks(-(sig - np.mean(sig)), height = thres)[0]
        p2 = find_peaks((sig - np.mean(sig)), height = thres)[0]
        if np.any(p1.tolist() + p2.tolist()):
            sig2[np.array(p1.tolist()+ p2.tolist())] = 1
    return sig2

def events(sample):
    """ Parameters
    --------
    sample : discretized array of recorded voltages
    
    Returns
    --------
    n : # of events in each temporal frame
    """
    if sample.ndim > 2:
        sample = sample.reshape(sample.shape[0],-1)
    n = np.array(np.sum(sample, axis = 1), dtype = int)
    return n 

def avinterv(n): 
    """ Parameters
    --------
    n : # of events in each temporal frame
    
    Returns
    --------
    avinterv : average inter event interval 
    """
    
    indici = np.arange(0,len(n),1)
    idx = indici[np.array(n,dtype = bool)]
    intertempi = idx[1:] - idx[:-1] 
    avinterv = int(round(np.mean(intertempi)))
    return avinterv


def binning(n, interv): 
    """
    Bins the data with the temporal bin "interv" and detects avalanches (one avalanche ends when an empty bin is found)
    
    Parameters
    n : # of events in each temporal frame (returned by the function "events")
    interv : temporal bin. Typical width is the ne of the average inter event interval (returned by "avinterv")
    --------
    Returns
    --------
    A list containing the events in avalanches. Each element of the returned list is an avalanche. The sum of the entries of the i-th avalanche is the size, the length of the i-th avalanche is its duration.
    """
    if len(n)%interv > 0:
    
        add = (int(len(n)/interv) + 1)* interv - len(n)
        n = n.tolist()
        for i in range(add):
            n = n + [0]
        
    n = np.asarray(n).reshape(int(len(n)/interv), interv)
    
    avalanches = []
    avalanches.append([])

    j = 0
    for z in range(len(n)):
        if np.any(n[z]):
            avalanches[j].append(np.sum(n[z]))
        else:
            j = j + 1
            avalanches.append([])
     
            
    for i in range(avalanches.count([])):
        avalanches.remove([])

    return avalanches

def draw(sample, xmin,xmax, model,ax= None,lim1 = 4,color = 'green'):
    """
    Fits the data contained in "sample" (i.e. sizes or durations of the avalanches) with the model chosen.
    xmin : the maximum xmin that can be considered in the power law fit.
    lim1 : upper limit of the range in which to search the power law parameter
    ax = object returned by fig.add_subplot()
    """
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
    ns = len(pwl.pdf(sample)[0])
    nbins = np.logspace(np.log10(min(sample)),np.log10(max(sample)),ns)
    nbins = nbins.tolist()
    ax.set_xscale('log')
    ypred = pwl.Fit(sample,xmin =(1,xmin + 1),xmax = xmax, parameter_range = {"alpha" : [1,lim1]},discrete = True)
    ax.hist(sample, density = True, histtype = 'bar',log = True,bins = nbins, color = color)
    pwl.plot_pdf(sample, color='r', linewidth=2, label='pdf', linear_bins = False)
    
    if model == 'power_law':
        ypred.power_law.plot_pdf( color='blue', linestyle='-', linewidth=2, label='Power Law fit')
        print(ypred.distribution_compare('power_law', 'exponential', normalized_ratio = True))
        print('Parameters are (alpha)',ypred.power_law.alpha)
    elif model == 'lognormal':
        ypred.lognormal.plot_pdf( color='blue', linestyle='-', linewidth=2, label='Lognormal fit')
        print('Parameters are (mu,sigma)',ypred.lognormal.mu,'and',ypred.lognormal.sigma)
    elif model == 'exponential':
        ypred.exponential.plot_pdf( color='blue', linestyle='-', linewidth=2, label='Exponential fit')
        print('Parameters are (lambda)',ypred.exponential.Lambda)
    elif model == 'stretched_exponential':
        ypred.stretched_exponential.plot_pdf( color='blue', linestyle='-', linewidth=2, label='Stretched_exponential fit')
        print('Parameters are (lambda,beta)',ypred.stretched_exponential.Lambda, 'and',ypred.stretched_exponential.beta)
    elif model == 'truncated_power_law':
        ypred.truncated_power_law.plot_pdf( color='blue', linestyle='-', linewidth=2, label='Truncated_power_law fit')
        print('Parameters are (alpha,lambda)',ypred.truncated_power_law.alpha, 'and',ypred.truncated_power_law.Lambda)
    plt.legend(fontsize = 'x-large')

def disegno(sizess, durationss, xmins, xmind, boool, maxs ,maxd,lim1,lim2):
   
    plt.figure()
    plt.xlabel("Avalanche size", fontsize = 'x-large')
    plt.ylabel("Probability density", fontsize = 'x-large')

    plt.xscale('log')
   
    ypred = pwl.Fit(sizess,xmin =(1,xmins+1),xmax = maxs, parameter_range = {'alpha' : [1,lim1]},discrete = True)
    
    new = []
    for i in range(len(sizess)):
        if sizess[i] >= ypred.xmin:
            new.append(sizess[i])
    
    ns = len(pwl.pdf(sizess)[0])
    

    nbins = np.logspace(np.log10(min(sizess)),np.log10(max(sizess)),ns)
    nbins = nbins.tolist()
   
    plt.hist(sizess, density = True, histtype = 'bar',log = True, bins = nbins,color = "green")
    pwl.plot_pdf(sizess, color='r', linewidth=2, label= 'pdf', linear_bins = False)
    ypred.power_law.plot_pdf( color='blue', linestyle='-', linewidth=2, label='Power Law fit')
    
    print(ypred.distribution_compare('power_law', 'exponential', normalized_ratio = True),
    ypred.distribution_compare('power_law', 'lognormal', normalized_ratio = True), ypred.power_law.xmin)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(fontsize = 'medium',  loc = 1)
 
    print('The predicted esponent and its error (sizes) are :', ypred.power_law.alpha, ' ', ypred.power_law.sigma)
    plt.figure()
    
    ypred2 = pwl.Fit(durationss,xmin =(1,xmind+1),xmax = maxd,parameter_range = {'alpha' : [1.,lim2]},discrete = boool)
    new2 = []
    for i in range(len(durationss)):
        if durationss[i] >= ypred2.xmin:
            new2.append(durationss[i])
            
    nd = len(pwl.pdf(durationss)[0])
    print(ns,nd)
 
    nbins2 = np.logspace(np.log10(min(durationss)),np.log10(max(durationss)),nd)
    nbins2 = nbins2.tolist()
    plt.xlabel("Avalanche duration", fontsize = 'x-large')
    plt.ylabel("Probability density", fontsize = 'x-large')
    plt.hist(durationss,density = True, log = True, bins = nbins2)
   
    pwl.plot_pdf(durationss, color='r', linewidth=2, label='pdf', linear_bins = False)
    ypred2.power_law.plot_pdf(color='blue', linestyle='-', linewidth=2, label='Power Law fit')
    print(ypred2.distribution_compare('power_law', 'exponential', normalized_ratio = True),
    ypred2.distribution_compare('power_law', 'lognormal', normalized_ratio = True), ypred2.power_law.xmin)
    plt.xscale('log')
    plt.legend(fontsize = 'medium', loc = 1)

    print('The predicted esponent and its error (durations) are :',ypred2.power_law.alpha, ' ', ypred2.power_law.sigma)
    
def RasterPlot(sample, av):
    """
    Parameters
    --------
    sample : Array of discretized data. Dimensions : temporal dim x spatial dim1 (x spatial dim2)
    av : width of temporal bin 
    
    Returns
    --------
    Plots the Raster Plots and the detected avalanches (an avalanche is preceded and followed by white bins)
    
    """
    
    sample = sample.reshape(sample.shape[0],-1) 
    s = binn(events(sample),av)
    times = np.asarray([av*i for i in range(len(s)+1)])
    s = np.asarray(s, dtype = bool)
  

    x= np.arange(0,sample.shape[1])
    for i in range(len(times)-1):
  
        if s[i]:
            plt.plot([times[i] for r in range(len(x))],x, color = 'red', linewidth = 1)
            plt.fill_betweenx(x,np.asarray([times[i] for r in range(len(x))]),np.asarray([times[i+1] for r in range(len(x))]), color = 'red', alpha = 0.3)

        else:
            plt.plot([times[i] for r in range(len(x))],x, '-', color = 'white', alpha = 0.3)
            plt.fill_betweenx(x,np.asarray([times[i] for r in range(len(x))]),np.asarray([times[i+1] for r in range(len(x))]), color = 'white')

    
    for i in range(sample.shape[0]):
        for j in range(sample.shape[1]):
            if sample[i,j]==1:
                plt.plot(i,j, 'b.')
    plt.xlabel('Time (temporal frame)', fontsize = 'xx-large')
    plt.ylabel('Channel',fontsize = 'xx-large')
    plt.xticks(fontsize = 'large')
    plt.yticks(fontsize = 'large')
    
def xminn(sample,xmin,xmax,lim):
    ypred = pwl.Fit(sample, xmin = (1,xmin+1),xmax = xmax, parameter_range = {'alpha' : [1.,lim]}, discrete = True)
    return ypred.power_law.xmin
    
def esponente(sample,xmin, boool,lim):
    ypred = pwl.Fit(sample,xmin = (1,xmin + 1),xmax = max(sample), parameter_range = {'alpha': [1,lim]}, discrete = boool)
    return ypred.power_law.alpha, ypred.power_law.sigma


def intertempi(n): 
    
    indici = np.arange(0,len(n),1)
    idx = indici[np.array(n,dtype = bool)]
  
    intertempi = idx[1:] - idx[:-1] 
    return intertempi

def scaling(sizess, durationss):
    tot = [durationss,sizess]

    numdur = np.sort(list(set(durationss)))
    sdit = [[] for i in range(len(numdur))]
    for i in range(len(numdur)):
        for z in range(len(durationss)):
            if tot[0][z] == numdur[i]:
                sdit[i].append(tot[1][z])
    medie = []
    std = []
    num = []
    for i in range(len(sdit)):
        medie.append(np.mean(sdit[i]))
        std.append(np.std(sdit[i])/np.sqrt(len(sdit[i])))
        num.append(len(sdit[i]))
        
    return numdur, medie,std, num

def prediz(alpha, tau):
    return (alpha - 1)/(tau - 1)

def ypred(sample, xmin, xmax,lim):
    """
    here xmin is a precise value, i. e. the best xmin previously selected by the fit
    """
    ypred = pwl.Fit(sample,xmin = (xmin,xmin + 1),xmax = xmax,parameter_range = {'alpha' : [1,lim]},discrete = True)
    return ypred

def goodness(sample, xmin, xmax, Ds):
 
    predic = ypred(sample, xmin, xmax)
    alpha = predic.power_law.alpha
    Dreal = predic.power_law.D
    y = np.arange(100)
    score = percentileofscore(Ds,Dreal)
    pval = 1-score/100.
    return pval

def modelcompare(sample,xmin,xmax,model1,model2):
    ypred = pwl.Fit(sample,xmin =(1,xmin+1),xmax = xmax,discrete = True)
    return ypred.distribution_compare(model1, model2, normalized_ratio = True)

def binningone(n):
    avalanches = []
    avalanches.append([])
    j = 0
    for z in range(len(n)):
        if n[z] > 0:
            avalanches[j].append(n[z])
        else:
            j = j + 1
            avalanches.append([])

    for i in range(avalanches.count([])):
        avalanches.remove([])
    return avalanches

def picchibinned(sample1, interv): 
    """
    sample1 : discretized data, dimensions: temporal dim x spatial dim1
    
    """
    if sample1.ndim > 2:
        sample1 = sample1.reshape(sample1.shape[0],-1)

    lenn = len(sample1)

    if lenn %interv != 0: 
    
        add = (int(len(sample1)/interv) + 1)* interv - len(sample1)
    
        addvec = np.zeros((add,len(sample1[0])))
        sample1 = np.vstack((sample1,addvec))
        
        for i in range(add):
                lenn = lenn + 1
            
    newvec = np.empty((int(lenn/interv),len(sample1[0])), dtype = int)
    prova = sample1.reshape(int(lenn/interv), interv,len(sample1[0]))

    
    for s in range(len(sample1[0])):
        for l in range(len(prova[:,:,s])):
            if np.any(prova[l,:,s]):
                newvec[l,s] = 1
            else:
                newvec[l,s] = 0
    
    n = np.zeros(len(newvec), dtype = int) 

    n[np.array(np.sum(newvec, axis = 1), dtype = bool)] = 1
    
    
    return n

def binn(n,interv):
    v = []
    add = (int(len(n)/interv) + 1)* interv - len(n)
    n = n.tolist()
    for i in range(add):
        n = n + [0]
        
    n = np.asarray(n).reshape(int(len(n)/interv), interv)
    
    for z in range(len(n)):
        if np.any(n[z]):
            v.append(1)
        else:
           
            v.append(0)
    return v
    

def Thres(coef,x):
    """
    Calculates Quiroga detection threshold
    """
    return coef*np.median(np.abs(x)/0.6745) 

    
def GaussianComparison(sample,n1,std,thres,n2 = 1):
    """
    Compares the distribution of the amplitudes in one electrode with the best Fit Gaussian
    
    Parameters
    --------
    sample : tridimensional array of the continuous data. Dimensions : temporal dim x spatial dim1 (x spatial dim2)
    (n1,n2) : coordinates of the chosen electrode
    """
    if sample.ndim < 3:
        sample = sample.reshape(sample.shape[0],sample.shape[1],1)
        std = std.reshape(sample1.shape[1],1)
    
    fig = plt.figure()

    
    plt.title('Electrode in coordinates (%d,%d) [Array dimensions: 4x55] ' %(n1+1,n2+1), fontsize  ='xx-large')

    sig = (sample[:,n1,n2] - np.mean(sample[:,n1,n2]))/std[n1,n2] # normalizes the signal by the Threshold (std/median)

    bins = np.arange(-6,6, 0.1)
    a,b = np.histogram(sig, bins =bins , density = True)
    x = (b[:-1]+ b[1:])/2 
    y = a
 
    n = len(x)                          
    means = sum(x*y)/n                   
    sigma = sum(y*(x-means)**2)/n        

    def gaus(x,x0,sig):
        return 1/(np.sqrt(2*np.pi*sig**2))*np.e**(-(x-x0)**2/(2*sig**2))

    popt,pcov = curve_fit(gaus,x,y,p0=[means,sigma])
 

    plt.hist(sig, bins = len(x), density = True)

    plt.plot(x,gaus(x,*popt),'-', color = 'purple',label='Best Gaussian fit')

    plt.xticks( np.arange(-6,6,1),fontsize = 'xx-large')
    plt.yticks( fontsize = 'xx-large')
    plt.xlabel('Normalized amplitude (SD)', fontsize = 'xx-large')
    plt.ylabel('Probability density', fontsize = 'xx-large')
    v = np.arange(0,0.6,0.01)
    x1 = np.array([-thres for i in range(len(v))])
    x2 = np.array([thres for i in range(len(v))])
    plt.plot(x1,v, 'r-', label = r'$\pm$ 3 SD')
    plt.plot(x2,v, 'r-', r'3 SD')
    plt.yscale('log')
    plt.legend( fontsize = 'xx-large')



    