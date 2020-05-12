import numpy as np
import matplotlib.pyplot as plt
import powerlaw as pwl
from matplotlib import cm
from statsmodels.regression import linear_model as sm
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import percentileofscore
from data_io_ts import *
from stats import *
from powerlaw_fit import *
# import setup for figure layout (see src/modules_plotting.py to customize)
from modules_plotting import * 
import matplotlib.cm as cm

def threshold(sample1,means,stds,thres,choose = "posneg", opz = "option1"):
    """ 
    Detects as events the points of maximum excursion over a threshold, considering either positive and negative excursions or only negative. if "option1" is selected, the one largest maximum between two crossings of the mean assigns the final event time.
    For a faster thresholding use the function below findpeaks.
    
    Parameters
    --------
    sample1 : tri or bidimensional array of recorded voltages (or even single time series from one electrode) (shape = temporal dim x spatial dim1 (x  spatial dim2) ) 
    means : array of the means of the signals (shape = spatial dim1 (x spatial dim2) ) 
    stds : array of the thresholds for each channel (standard deviations/medians...) (shape = spatial dim1 (x spatial dim2)
    thres : multiplicative coefficent for the thresholds
    #ref : refractory period. Minimum distance between consecutive events.
    choose : if "posneg" both positive and negative deflections are considered, if "neg" only negative
    opz : if "option1" the one largest maximum between two crossings of the mean assigns the final event time, if "option2" an event is imply the point of maximum excursion over a threshold
    
    Returns
    --------
    sample2 : discretized array with the initial shape 
    """
    initshape = sample1.shape
    
    if sample1.ndim > 2:
        sample1 = sample1.reshape(sample1.shape[0],-1)
        means = means.reshape(-1,)
        stds = stds.reshape(-1,)
        
        
    if sample1.ndim == 1: # so this same code works even when considering a single time series
        sample1 = sample1.reshape(sample1.shape[0],1,1)
        means = means.reshape(1,1)
        stds = stds.reshape(1,1)
        
  
    
    if opz == "option1":
        sample2 = np.zeros(sample1.shape, dtype = int)
        if choose == "posneg":

            for s in range(sample1.shape[1]):
                if stds[s]>0:
                    sig = sample1[:,s]

                    stan = stds[s]
                    tempi = np.arange(0,len(sig),1)
                    prova = np.where(np.abs(sig - means[s]) >= thres*stan,True,False)
                    changesign1 = np.where(np.sign(sig[:-1]-means[s]) !=np.sign(sig[1:]-means[s]), True,False).reshape(-1,)
                    changesign = np.hstack((changesign1,False))

                    init = []
                    end = []

                    termined = True
                    started = False

                    if prova[0] == True:
                        init.append(0)
                        termined = False
                        started = True

                    for i in range(1,len(prova)):

                        if termined:

                            if prova[i-1] == False and prova[i] == True:
                                init.append(i)
                                termined = False
                                started = True


                        if changesign[i] and started:

                            end.append(i + 1)
                            termined = True
                            started = False




                    if termined == False:
                        end.append(len(prova)+1)


                    groups = []
                    times = []
                    for l in range(len(init)):
                        groups.append(np.abs(sig[init[l]:end[l]]-means[s]))

                        times.append(tempi[init[l]:end[l]])

                    zeta = []
                    for m in range(len(groups)):
                        zeta.append(times[m][groups[m].tolist().index(max(groups[m]))])

                    sample2[zeta,s] = 1


        if choose == "neg":

            for s in range(sample1.shape[1]):
                 if stds[s]>0:
                    sig = sample1[:,s]

                    stan = stds[s]
                    tempi = np.arange(0,len(sig),1)
                    prova = np.where((sig - means[s]) <= -thres*stan,True,False)
                    changesign1 = np.where(np.sign(sig[:-1]-means[s]) !=np.sign(sig[1:]-means[s]), True,False).reshape(-1,)
                    changesign = np.hstack((changesign1,False))

                    init = []
                    end = []

                    termined = True
                    started = False

                    if prova[0] == True:
                        init.append(0)
                        termined = False
                        started = True

                    for i in range(1,len(prova)):

                        if termined:

                            if prova[i-1] == False and prova[i] == True:
                                init.append(i)
                                termined = False
                                started = True


                        if changesign[i] and started:

                            end.append(i + 1)
                            termined = True
                            started = False




                    if termined == False:
                        end.append(len(prova)+1)


                    groups = []
                    times = []
                    for l in range(len(init)):
                        groups.append(np.abs(sig[init[l]:end[l]]-means[s]))

                        times.append(tempi[init[l]:end[l]])

                    zeta = []
                    for m in range(len(groups)):
                        zeta.append(times[m][groups[m].tolist().index(max(groups[m]))])

                    sample2[zeta,s] = 1


    if opz == "option2":
        sample2 = np.zeros(sample1.shape, dtype = int)
        if choose == "posneg":

            for s in range(sample1.shape[1]):
                if stds[s]>0:
                    sig = sample1[:,s]

                    stan = stds[s]
                    tempi = np.arange(0,len(sig),1)
                    prova = np.where(np.abs(sig - means[s]) >= thres*stan,True,False)

                    init = []
                    end = []
                    termined = True
                    started = False
                    if prova[0] == True:
                        init.append(0)
                        termined = False
                        started = True
                    for i in range(1,len(prova)):

                        if prova[i-1] == False and prova[i] == True and not started:
                            init.append(i)
                            termined = False
                            started = True
                        if (i+1) < len(prova):
                            if prova[i] == True and prova[i+1] == False and started:

                                end.append(i + 1)
                                termined = True
                                started = False 



                    if termined == False:
                        end.append(len(prova)+1)


                    groups = []
                    times = []
                    for l in range(len(init)):
                        groups.append(np.abs(sig[init[l]:end[l]]-means[s]))

                        times.append(tempi[init[l]:end[l]])

                    zeta = []
                    for m in range(len(groups)):
                        zeta.append(times[m][groups[m].tolist().index(max(groups[m]))])

                    sample2[zeta,s] = 1


        if choose == "neg":

            for s in range(sample1.shape[1]):
                if stds[s]>0:
                    sig = sample1[:,s]

                    stan = stds[s]
                    tempi = np.arange(0,len(sig),1)
                    prova = np.where((sig - means[s]) <= -thres*stan,True,False)

                    init = []
                    end = []
                    termined = True
                    started = False

                    if prova[0] == True:
                        init.append(0)
                        termined = False
                        started = True
                        #print('started',i)

                    for i in range(1,len(prova)):

                        if prova[i-1] == False and prova[i] == True and not started:
                            init.append(i)
                            termined = False
                            started = True
                            #print('started', i)
            
                            
                        if (i+1) < len(prova):
                            if prova[i] == True and prova[i+1] == False and started:
                                end.append(i + 1)
                                termined = True
                                started = False
                                #print('end',i)

                    if termined == False:
                        end.append(len(prova)+1)
                        #print('end', len(prova)+1)
                        
                    #print(len(init) == len(end))
                    #print(len(init))
                    #print(len(end))
                    groups = []
                    times = []
                    for l in range(len(init)):
                        groups.append(np.abs(sig[init[l]:end[l]]-means[s]))

                        times.append(tempi[init[l]:end[l]])

                    zeta = []
                    for m in range(len(groups)):
                        zeta.append(times[m][groups[m].tolist().index(max(groups[m]))])

                    sample2[zeta,s] = 1

    return sample2.reshape(initshape)


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
    idx = np.where(n > 0)[0]
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
    
    new = np.array(np.sum(n,axis = 1),dtype = int)
    
    init = []
    end = []
 
    prova = new>0
    
    if prova[0] == True:
        init.append(0)
        
    
        
    for i in range(1,len(prova)):
    
        if prova[i-1] == False and prova[i] == True:
            init.append(i)
        if prova[i-1] == True and prova[i] == False:
            end.append(i)
    
    if prova[-1] == True:
        end.append(len(prova)+1)
            
    avalanches = []

    for s in range(len(init)):
        avalanches.append(new[init[s]:end[s]])
    
    sizes = []
    durations = []
    for l in range(len(avalanches)):
        sizes.append(np.sum(avalanches[l]))
        durations.append(len(avalanches[l]))
        
    return sizes, durations
    
    
def avalanches(n, interv): 
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
    
    new = np.array(np.sum(n,axis = 1),dtype = int)
    tempi = np.arange(0,len(new),1)
    times = []
    init = []
    end = []
 
    prova = new>0
    
    if prova[0] == True:
        init.append(0)
        
    
        
    for i in range(1,len(prova)):
    
        if prova[i-1] == False and prova[i] == True:
            init.append(i)
        if prova[i-1] == True and prova[i] == False:
            end.append(i)
    
    if prova[-1] == True:
        end.append(len(prova)+1)
            
    avalanches = []

    for s in range(len(init)):
        avalanches.append(new[init[s]:end[s]])
        times.append(tempi[init[s]:end[s]])
    
    sizes = []
    durations = []
    for l in range(len(avalanches)):
        sizes.append(np.sum(avalanches[l]))
        durations.append(len(avalanches[l]))
        
    return sizes, durations, avalanches, times
    
    
"""
def binning(n,interv)
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
"""
def main(sample,av = "default"):
    """
    Parameters
    --------
    sample : discretized array of voltages
    av : width of the interval to bin the data and calculate avalanches (average inter event interval usually)
    
    Returns
    --------
    sizes and durations of the detected avalanches
    """
    peaks = events(sample)
    if av == "default":
        av = avinterv(peaks)
    sizes, durations = binning(peaks,av)
    
    return sizes, durations

def draw(sample, xmin,xmax, model,ax= None,lim1 = 4,color = 'green',normalized = False):
    """
    Fits the data contained in "sample" (i.e. sizes or durations of the avalanches) with the model chosen.
    xmin : the maximum xmin that can be considered in the power law fit.
    lim1 : upper limit of the range in which to search the power law parameter
    ax = object of the type AxesSubplot (returned by fig.add_subplot())
    """
    ypred = pwl.Fit(sample,xmin =(1,xmin + 1),xmax = xmax, parameter_range = {"alpha" : [1,lim1]},discrete = True)
    
    if normalized == True:
        new = []
        for i in range(len(sample)):
            if sample[i] >= ypred.xmin:
                new.append(sample[i])
                
        sample = new
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
    ns = len(pwl.pdf(sample)[0]) -1# == nbins + 1
    print(ns)
    
    ax.set_xscale('log')
    

    nbins = np.logspace(np.log10(min(sample)),np.log10(max(sample)),ns)
    ax.hist(sample, density = True, histtype = 'stepfilled',log = True,bins = nbins, color = color, alpha = 0.25)
    pwl.plot_pdf(sample, color='r', linewidth=2, label='pdf', linear_bins = False)
    
    if model == 'power_law':
        ypred.power_law.plot_pdf( color='blue', linestyle='-', linewidth=2, label='Power Law fit')
        print(ypred.distribution_compare('power_law', 'exponential', normalized_ratio = True))
        print('Parameters are (alpha)',ypred.power_law.alpha, 'xmin',ypred.power_law.xmin)
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
    plt.legend()

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
    
    ns = len(pwl.pdf(sizess)[0])-1
    

    nbins = np.logspace(np.log10(min(sizess)),np.log10(max(sizess)),ns)
    
   
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
            
    nd = len(pwl.pdf(durationss)[0])-1
    print(ns,nd)
 
    nbins2 = np.logspace(np.log10(min(durationss)),np.log10(max(durationss)),nd)
    
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

    
    for j in range(sample.shape[1]):
        idx = np.where(sample[:,j]> 0)[0]
        plt.plot(idx,[j for i in range(len(idx))], 'b.')
        
    plt.xlabel('Time (temporal frame)')
    plt.ylabel('Channel')
    plt.xticks()
    plt.yticks()
    


def intertempi(n): 
    idx = np.where(n > 0)[0]
    intertempi = idx[1:] - idx[:-1] 
    return intertempi

def sgivent(sizess, durationss):
    tot = [durationss,sizess]

    numdur = np.sort(list(set(durationss)))
    sdit = [[] for i in range(len(numdur))]
    for i in range(len(numdur)):
        for z in range(len(durationss)):
            if tot[0][z] == numdur[i]:
                sdit[i].append(tot[1][z])
    medie = []
    std = []

    for i in range(len(sdit)):
        medie.append(np.mean(sdit[i]))
        std.append(np.std(sdit[i])/np.sqrt(len(sdit[i])))
        
    return numdur,medie,std

def delta(alpha, salpha, tau, stau):
    return (alpha - 1)/(tau -1), np.sqrt((1/(tau -1))**2*salpha**2+ ((1-alpha)/(tau -1)**2)**2*stau**2)


def scaling(sizes, durations, avinterval, ax = None, tau = "default", errtau = "default", alpha = "default", erralpha = "default", maxxminsizes = "default", maxxmindur = "default", xmaxsizes = "default", xmaxdur = "default", lim1 = 4 , lim2 = 4):
    
    
    if ax == None:
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot(1,1,1)
        
    if maxxminsizes ==  "default":
        maxxminsizes = max(sizes)
    if xmaxsizes ==  "default":
        xmax = max(sizes)

    
    if maxxmindur ==  "default":
        maxxmindur = max(durations)
    if xmaxdur ==  "default":
        xmaxdur = max(durations)
    
    
    if tau == "default" and errtau == "default":
        tau,errtau = esponente(sizes,maxxmin = maxxminsizes,xmax =xmaxsizes,lim = lim1)
        
    if alpha == "default" and erralpha == "default":
        alpha,erralpha = esponente(durations,maxxmin = maxxmindur,xmax =xmaxdur,lim = lim2)

    
    pred =  delta(alpha, erralpha, tau, errtau )[0]
    errpred = delta(alpha, erralpha,tau, errtau)[1]


    durations = np.array(durations)*avinterval
    
    xmin1 = min(sizes)
    xmin2 = min(durations)
    
    
    prova = np.array([np.asarray(sizes), np.asarray(durations)])
    prova = prova.transpose()
    prova2 = [0 for i in range(len(prova))]
    
    for r in range(len(prova)):
        if prova[r][0] < xmin1 or prova[r][1] < xmin2 :
            prova2[r] = False
        else:
            prova2[r] = True
            
    new = prova[prova2]
    
    a,b,c = sgivent(new[:,0],new[:,1])
    x = np.hstack((np.log10(a).reshape(-1,1), np.ones(len(a)).reshape(-1,1)))
    
    y = np.log10(b).reshape(-1,1)

    ols = sm.OLS(y,x)

    ols_result = ols.fit()
    
    fit = ols_result.params[0]
    errfit = ols_result.bse[0]
    inter =ols_result.params[1]

    grays = cm.Greys(np.linspace(0,1,15))
    reds = cm.Reds(np.linspace(0,1,15))
    greens = cm.Greens(np.linspace(0,1,15))
    blues = cm.Blues(np.linspace(0,1,30))

    nott = [0 for i in range(len(prova))]
    for r in range(len(prova)):
        if prova[r][0] < xmin1 or prova[r][1] < xmin2 :

            nott[r] = True
        else:

            nott[r] = False

    print('Prediction from crackling noise relation: delta = ',pred, '+-', errpred)
    print('Fit from of average size given duration points: delta = ',fit, '+-', errfit)
    
 

    ax.set_xscale('log')
    ax.set_yscale('log')
   
    x = np.arange(min(durations),max(durations))

    ax.plot(durations, sizes, '.', color = blues[15],alpha = 0.3)
    ax.plot(np.asarray(durations)[nott],np.asarray(sizes)[nott], '.', color = 'gray')
    
    ax.plot(x, (10**inter)*x**pred, 'r', label = 'Prediction', lw = 1.5)

    ax.plot(x, (10**inter)*x**fit, color = greens[12], label = 'Fit', lw = 1.5)

 
    if errfit > errpred:
        ax.fill_between(x,(10**(errfit*-3))*(10**inter)*x**fit,
                         (10**(errfit*3))*(10**inter)*x**fit,color= greens[5],lw=0)


        ax.fill_between(x,(10**(errpred*-3))*(10**inter)*x**pred,
                         (10**(errpred*+3))*(10**inter)*x**pred,color= reds[5],lw=0)
    
    else:
        ax.fill_between(x,(10**(errpred*-3))*(10**inter)*x**pred,
                         (10**(errpred*+3))*(10**inter)*x**pred,color= reds[5],lw=0)
        
        ax.fill_between(x,(10**(errfit*-3))*(10**inter)*x**fit,
                         (10**(errfit*3))*(10**inter)*x**fit,color= greens[5],lw=0)


    
        
    ax.set_xlabel('Durations [ms]')
    ax.set_ylabel('Sizes [# of events]')


    ax.legend()

    ax.errorbar(a,b,yerr = c, fmt = 'o', color = blues[29],markersize = 3.5,barsabove = False,capsize = 3, elinewidth = 3,capthick = 1)
    

def xminn(sample,maxxmin = "default",xmax = "default",lim = 4):
    if maxxmin ==  "default":
        maxxmin = max(sample)
    if xmax ==  "default":
        xmax = max(sample)
  
    ypred = pwl.Fit(sample, xmin = (1,maxxmin+1),xmax = xmax, parameter_range = {'alpha' : [1.,lim]}, discrete = True)
    return ypred.power_law.xmin
    

def esponente(sample,maxxmin = "default",xmax = "default",lim = "default"):
    if maxxmin ==  "default":
        maxxmin = max(sample)
    if xmax ==  "default":
        xmax = max(sample)
    if lim ==  "default":
        lim = 4
        
    ypred = pwl.Fit(sample,xmin = (1,maxxmin + 1),xmax = xmax, parameter_range = {'alpha': [1,lim]}, discrete = True)
    return ypred.power_law.alpha, ypred.power_law.sigma

def ypred(sample, xmin, xmax,lim = 4):
    """

    """
    #xmin qui è un valore preciso
    ypred = pwl.Fit(sample,xmin = (xmin,xmin + 1),xmax = xmax,parameter_range = {'alpha' : [1,lim]},discrete = True) # attenzione che qui ho modificato
    return ypred

def goodness(sample, xmin, xmax, Ds,lim):
    #mi ricordo che lo ho modificato
 
    predic = ypred(sample, xmin, xmax,lim)
    alpha = predic.power_law.alpha
    print(alpha)
    Dreal = predic.power_law.D
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
    prova = np.sum(sample1.reshape(int(lenn/interv), interv,len(sample1[0])),axis = 1)
    #verificare
    
    for s in range(len(sample1[0])):
        for l in range(len(prova)):
            if (prova[l,s]):
                newvec[l,s] = 1
            else:
                newvec[l,s] = 0
    
    n = np.zeros((len(newvec)), dtype = int) 

    n[np.array(np.sum(newvec, axis = 1), dtype = bool)] = 1
    
    
    return n

def binn(n,interv):
    v = []
    
    if len(n)%interv > 0:
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

    
def GaussianComparison(sample,Thres,coef, n1, n2 = 0):
    """
    Compares the distribution of the amplitudes in one electrode with the best Fit Gaussian
    
    Parameters
    --------
    sample : tri/bidimensional array of the continuous data. Dimensions : temporal dim x spatial dim1 (x spatial dim2)
    Thres : tri/bidimensional array of the standard deviations of the data/medians of the data... Dimensions : temporal dim x spatial dim1 (x spatial dim2)
    coef : multiplicative coefficient for the thresholds
    (n1,n2) : coordinates of the chosen electrode
    """
    if sample.ndim < 3:
        sample = sample.reshape(sample.shape[0],sample.shape[1],1)
    if Thres.ndim < 3:
        Thres = Thres.reshape(Thres.shape[0],1)
        
    fig = plt.figure()

    
    plt.title('Electrode in coordinates (%d,%d) [Array dimensions: 4x55] ' %(n1+1,n2+1))

    sig = (sample[:,n1,n2] - np.mean(sample[:,n1,n2]))/(Thres[n1,n2]) # normalizes the signal by the SD

    bins = np.arange(min(sig),max(sig)+0.1, 0.1)
    a,b = np.histogram(sig, bins =bins , density = True)
    x = (b[:-1]+ b[1:])/2 
    y = a

    popt,pcov = curve_fit(gaus,x,y,p0=[0,1],bounds = (np.array([0-0.5, -4]),np.array([0+0.5,4])))
 
    plt.hist(sig, bins = bins, density = True)

    plt.plot(x,gaus(x,*popt),'-', color = 'purple',label='Best Gaussian fit')

    plt.xticks( np.arange(int(min(sig)),int(max(sig)),1))
    plt.yticks( fontsize = 'large')
    plt.xlabel('Normalized amplitude (SD)')
    plt.ylabel('Probability density')
    v = np.arange(0,max(y),0.01)
    x1 = np.array([-coef for i in range(len(v))])
    x2 = np.array([coef for i in range(len(v))])
    plt.plot(x1,v, 'r-', label = r'$\pm$ %s SD' %coef)
    plt.plot(x2,v, 'r-', r'%s SD' %coef)
    plt.yscale('log')
    plt.legend()
    plt.ylim(10**-7,)

    
def GaussianComparisonMean(sample, xlim1,xlim2,thres):

    pdf = []

    if sample.ndim > 2:
        sample = sample.reshape(-1,sample.shape[1]*sample.shape[2])

    for s in range(sample.shape[1]):
        if np.std(sample[:,s]) > 0:
            sig = (sample[:,s]- np.mean(sample[:,s]))/np.std(sample[:,s])
            bins = np.arange(xlim1,xlim2+0.1, 0.1)
            a,b = np.histogram(sig, bins =bins , density = True)
            pdf.append(a)
    coord = (b[:-1]+ b[1:])/2 

    for s in range(219):
        plt.plot(coord, pdf[s], '.-',color = 'gray', alpha = 0.1)
    pdfmean = np.mean(pdf,axis = 0)
    pdfstd = np.std(pdf,axis = 0)/np.sqrt(len(pdf))

    x = coord
    y = pdfmean


    popt,pcov = curve_fit(gaus,x,y,p0=[0,1],bounds = (np.array([0-0.5, -4]),np.array([0+0.5,4])))



    plt.plot(x,gaus(x,*popt),'-', color = 'purple',label='Best Gaussian Fit')


    plt.xlabel('Normalized Amplitudes (SD)')
    plt.ylabel('Probability Density')
    v = np.arange(0,0.6,0.01)

    x1 = np.array([-thres for i in range(len(v))])
    x2 = np.array([thres for i in range(len(v))])
    plt.plot(x1,v, 'r-', label = r'$\pm$ %s SD' %thres)
    plt.plot(x2,v, 'r-', r'%s SD' %thres)
    plt.errorbar(x,y,yerr = pdfstd, fmt = '-', color = 'blue',ecolor = 'blue', elinewidth = 1,capsize = 3, barsabove = True, label = 'Average')
    plt.plot(x,y, 'b-')
    plt.yscale('log')
    plt.xticks(np.arange(xlim1,xlim2+1,1))
    plt.legend()


def LogScript(x_data,av,x_datarnd,avrnd,legend = True, ax = None,color = 'blue',nbins = "default", nbinsrnd = "default"):   
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
    xmin = xminn(x_data,max(x_data),max(x_data),4)
    xmax = max(x_data)
    
    new1 = []

    for i in range(len(x_data)):
        if x_data[i] >= xmin:
            new1.append(x_data[i])
        
    if nbins == "default":
        nbins = len(pwl.pdf(new1)[0])-1
        
    print('NBins for non random data are',len(pwl.pdf(new1)[0])-1)
    
    nrep_synth = 0
    x,nx = xdata_to_xnx(new1,norm=False,xmin=xmin,xmax=xmax)
    result = fit_power_disc_sign(x,nx,xmin=xmin,xmax=xmax,nrep_synth=nrep_synth)
    alpha = result['alpha']
    new = np.array(new1)
    new = new*av
    xmin = xmin*av
    xmax = xmax*av
    x,nx = xdata_to_xnx(new,norm=False,xmin=xmin,xmax=xmax)
    px_fit = pdf_power_disc(x,xmin,xmax,alpha) 
    
  
        
    
    bins = 10**(np.arange(min(np.log10(new1)), max(np.log10(new1)) +(max(np.log10(new1))-min(np.log10(new1)))/nbins,(max(np.log10(new1))-min(np.log10(new1)))/nbins))
    
    bins = bins*av
    
    pdf = [[] for i in range(len(bins)-1)]
    for i in range(1,len(bins)):
        for l in range(len(new)):
            if new[l]>= bins[i-1] and new[l] < bins[i]:
                pdf[i-1].append(new[l])

    pdfnorm = [len(pdf[i-1])/((bins[i]-bins[i-1])*len(new))for i in range(1,len(bins))]


    centredbin = [(bins[i]+ bins[i + 1])/2 for i in range(len(bins)-1)]

#------------------Random Data----------------------------------------------------------------------------------#

    datarnd = np.array(x_datarnd)*avrnd
    xmin2 =xmin
    new2 = []

    for i in range(len(datarnd)):
        if datarnd[i] >= xmin2:
            new2.append(datarnd[i])

    new2 = np.array(new2) 
    if nbinsrnd == "default":
        nbinsrnd = len(pwl.pdf(new2/avrnd)[0])-1
    print(len(pwl.pdf(new2/avrnd)[0])-1)

    bins = 10**(np.arange(min(np.log10(new2)), max(np.log10(new2)) +(max(np.log10(new2))-min(np.log10(new2)))/nbinsrnd,(max(np.log10(new2))-min(np.log10(new2)))/nbinsrnd))
    pdf = [[] for i in range(len(bins)-1)]
    for i in range(1,len(bins)):
        for l in range(len(new2)):
            if new2[l]>= bins[i-1] and new2[l] < bins[i]:
                pdf[i-1].append(new2[l])

    pdfnormrnd = [len(pdf[i-1])/((bins[i]-bins[i-1])*len(new2))for i in range(1,len(bins))]


    centredbinrnd = [(bins[i]+ bins[i + 1])/2 for i in range(len(bins)-1)]



    ax.plot(centredbin, pdfnorm,'ko-',  label = 'Data')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(x, px_fit, '-',color = color,label = 'Fit')
  

    ax.plot(centredbinrnd, pdfnormrnd,'--', color = 'gray', label = 'Shuffled')

#pwl.plot_pdf(new, color = 'black', marker = 'o')
    if legend:
        ax.legend(loc = (1,0),framealpha = 0)
    #plt.savefig('Durationsspikesrest.png', dpi = 300, bbox_inches ='tight', transparent = True)
    
def RocCurve(x, good_pdf, bad_pdf, coef = 2,coef2 = 0,reverse = False):
    multiplicative = True
    if reverse == True:
        bad_pdf = bad_pdf[::-1]
        good_pdf = good_pdf[::-1]
        x = x[::-1]
        
    #Total Bad
    fp = []
    for l in range(len(x)):
        if good_pdf[l]> 0 and good_pdf[l]<= coef*bad_pdf[l]:
            fp.append(good_pdf[l])
      
    total_bad = np.sum(fp)

    #Total Good (nota il fattore moltiplicativo coef)
    tp = []
    for l in range(len(x)):
        if multiplicative:
            if good_pdf[l]>coef*bad_pdf[l]:
                tp.append(good_pdf[l])
                
        else:
            if good_pdf[l]>coef2+bad_pdf[l]:
                tp.append(good_pdf[l])
            
    total_good = np.sum(tp)
   
    TPR_list=[]
    FPR_list=[]
      
    for i in x: #ciclo sulle possibili soglie
        
        if total_bad> 0 and total_good >0:
            
            #False Positive
            summ = []
            for l in range(len(x[:x.tolist().index(i)+1])):
                if good_pdf[l]<= coef*bad_pdf[l] and good_pdf[l]> 0:
                    summ.append(good_pdf[l])
                    
            FPR_list.append(np.sum(summ)/total_bad)
                    
            #True Positive
            summ2 = []
    
            for l in range(len(x[:x.tolist().index(i)+1])):
                if multiplicative:
                    if good_pdf[l]>coef*bad_pdf[l]:
                        summ2.append(good_pdf[l])
                
                else:
                    if good_pdf[l]>coef2+bad_pdf[l]:
                        summ2.append(good_pdf[l])
                    
            TPR_list.append(np.sum(summ2)/total_good)
            
            
    return FPR_list,TPR_list

def gaus(x,x0,sigm):
    return 1/(np.sqrt(2*np.pi*sigm**2))*np.e**(-(x-x0)**2/(2*sigm**2))

def OptimalThreshold(a1, coef):
    if a1.ndim > 2:
        a1 = a1.reshape(a1.shape[0],-1) #a1 matrice dei segnali. Ora ha shape temp dim x spatila dim (220)
        
    posthres = []
    negthres = []
    for s in range(a1.shape[1]): #ciclo sugli elettrodi, ora mi concentro su elettrodo 150 per farti vedere
        if np.std(a1[:,s])> 0: #elettrodo 15 non registra attività e non è da considerare

            sig = (a1[:,s] - np.mean(a1[:,s]))/np.std(a1[:,s]) # normalizes the signal by the SD

            bins = np.arange(min(sig),max(sig) + 0.1, 0.1)
            good_pdf,b = np.histogram(sig, bins =bins , density = True)

            x = (b[:-1]+ b[1:])/2 

            popt,pcov = curve_fit(gaus,x,good_pdf,p0=[0,1], bounds = (np.array([-0.5,-4]), np.array([0.5,4])))
            bad_pdf = gaus(x,*popt)


            idx = np.where(x> 0.00000001)[0][0] #MI concentro solo sulla parte sinistra delle distribuzioni (es:[-10,0])

            coord = x[idx:]
            coef2 = 0
 
            l, j = RocCurve(coord, good_pdf[idx:],bad_pdf[idx:],coef,coef2, True)
            #Mi ritorna i False Positive e i True negative rate, al variare di tutte le possibili soglie (coord)
            #Mi calcolo le distanze dalla diagonale, dato che il punto massimamente distante dalla diagonale è quello
            #ottimale


            dist = []
            new = []
            jnew = []
            lnew =[]
            for i in range(len(l)):
                if  j[i] >= l[i]:
                    dist.append(np.sqrt((l[i]-l[i])**2 +(j[i]-l[i])**2)) 
                    new.append(coord[::-1][i])
                    lnew.append(l[i])
                    jnew.append(j[i])

            if len(dist)> 0:
                posthres.append(new[dist.index(max(dist))]) #soglia corrispondene alla massima distanza
                
            if len(posthres) < 1:
                posthres.append(max(sig))
                


            #-------------negthres---------#

            sig = (a1[:,s] - np.mean(a1[:,s]))/np.std(a1[:,s]) # normalizes the signal by the SD

            bins = np.arange(min(sig),max(sig)+0.1, 0.1)
            good_pdf,b = np.histogram(sig, bins =bins , density = True)

            x = (b[:-1]+ b[1:])/2 


         

            popt,pcov = curve_fit(gaus,x,good_pdf,p0=[0,1], bounds = (np.array([-0.5,-4]), np.array([0.5,4])))
            bad_pdf = gaus(x,*popt)


            idx = np.where(x> 0.00000001)[0][0] #MI concentro solo sulla parte sinistra delle distribuzioni (es:[-10,0])


            coord = x[:idx]
            #0.0012
             #coef = 9.999999999e-05
            #coef = 1.5
            l, j = RocCurve(coord, good_pdf[:idx],bad_pdf[:idx],coef,coef2, False)
            #Mi ritorna i False Positive e i True negative rate, al variare di tutte le possibili soglie (coord)
            #Mi calcolo le distanze dalla diagonale, dato che il punto massimamente distante dalla diagonale è quello
            #ottimale


            dist = []
            new = []
            jnew = []
            lnew =[]
            for i in range(len(l)):
                if  j[i] >= l[i]:
                    dist.append(np.sqrt((l[i]-l[i])**2 +(j[i]-l[i])**2)) 
                    new.append(coord[i])
                    lnew.append(l[i])
                    jnew.append(j[i])

            if len(dist)> 0:
                negthres.append(new[dist.index(max(dist))]) #soglia corrispondene alla massima distanza
                
            if len(negthres) < 1:
                negthres.append(min(sig))
                
   
    print('Optimal positive threshold is ', np.mean(posthres)) #media sugli elettrodi
        #quando non me la ritorna perdo statistica?

      
    print('Optimal negative threshold is ', np.mean(negthres)) #media sugli elettrodi
 
        
    print('Optimal threshold is',np.mean(np.abs(np.array(negthres)).tolist()+ posthres), 'SD')
    return posthres,negthres,l,j,lnew,jnew,dist.index(max(dist))
#np.mean(np.abs(np.array(negthres)).tolist()+ posthres)



def DetectedEvents(a1,reverse,coef, s):
    #Visualizzo quali punti mi ha identificato come True Positive in elettrodo identificato da s
    """
    s : chosen electrode
    """
    multiplicative = True
    posthres,negthres,fpos,tpos,fp2,tp2,ind = OptimalThreshold(a1[:,s:s+1],coef)
    if reverse and len(posthres) > 0:
        i = posthres[0]
    if not reverse:
        i = negthres[0]

    
    sig = (a1[:,s] -np.mean(a1[:,s]))/np.std(a1[:,s]) # normalizes the signal by the SD

    bins = np.arange(min(sig),max(sig)+0.1,0.1)
    y,b = np.histogram(sig, bins =bins , density = True)
    x = (b[:-1]+ b[1:])/2 

    good_pdf = y

    popt,pcov = curve_fit(gaus,x,y,p0=[0,1], bounds = (np.array([-0.5,-4]), np.array([0.5,4])))
    bad_pdf = gaus(x,*popt)


    idx = np.where(x> 0.00000001)[0][0]

    if reverse == True:
        coord = x[idx:]
        x = coord        
        good_pdf = good_pdf[idx:]
        bad_pdf =  bad_pdf[idx:]

    else: 
        coord = x[:idx]
        x = coord        
        good_pdf = good_pdf[:idx]
        bad_pdf =  bad_pdf[:idx]

    if reverse == True:
        bad_pdf = bad_pdf[::-1]
        good_pdf = good_pdf[::-1]
        x = x[::-1]


    fp = []
    for l in range(len(x)):
        if good_pdf[l]> 0 and good_pdf[l]<= coef*bad_pdf[l]:
            fp.append(good_pdf[l])

    total_bad = np.sum(fp)
    FPR_list =[]
    TPR_list = []
    tp = []
    for l in range(len(x)):
        if multiplicative:
            if good_pdf[l]>coef*bad_pdf[l]:
                tp.append(good_pdf[l])
        else:
            if good_pdf[l]>=coef2+bad_pdf[l]:
                tp.append(good_pdf[l])

    good_x = []
    bad_x = []
    summ2 = []
    summ = []
    total_good = np.sum(tp)
    if total_bad> 0 and total_good >0:

            #false positive
        for l in range(len(x[:x.tolist().index(i)+1])):
            if good_pdf[l]<= coef*bad_pdf[l] and good_pdf[l]> 0:
                summ.append(good_pdf[l])
                bad_x.append(x[l])

        FPR_list.append(np.sum(summ)/total_bad)



        #True positive

        for l in range(len(x[:x.tolist().index(i)+1])):
            if multiplicative == True:
                if good_pdf[l]>coef*bad_pdf[l]:
                    summ2.append(good_pdf[l])
                    good_x.append(x[l])
            else:
                if good_pdf[l]>=coef2+bad_pdf[l]:
                    summ2.append(good_pdf[l])
                    good_x.append(x[l])

        TPR_list.append(np.sum(summ2)/total_good)
        
    plt.figure()
    plt.title('Roc Curve')
    plt.plot(fpos,tpos,'r-')
    plt.plot(fp2[ind],tp2[ind], 'ko', markersize = 10)
    plt.plot(fpos,fpos, 'b-')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    
    plt.figure()
    if reverse== False:
        plt.plot(coord, good_pdf, 'b-')
        plt.plot(good_x,summ2, 'ro', label = 'Identified positive')
        plt.plot(coord, bad_pdf, '.')

        plt.yscale('log')
        plt.legend()



    if reverse == True:
        plt.plot(x[::-1], good_pdf[::-1], 'b-')
        plt.plot(good_x,summ2,'ro', label = 'Identified positive')
        plt.plot(x[::-1], bad_pdf[::-1], '.')

        plt.yscale('log')
        plt.legend()

