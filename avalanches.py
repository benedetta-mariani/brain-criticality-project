import numpy as np
import matplotlib.pyplot as plt
import powerlaw as pwl
from matplotlib import cm
from statsmodels.regression import linear_model as sm
from scipy import signal 
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import percentileofscore
import cmath
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
                        end.append(len(prova))


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
                        end.append(len(prova))


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
                        end.append(len(prova))##


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
                        end.append(len(prova))
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
    interv : temporal bin. Typical width is the one of the average inter event interval (returned by "avinterv")
    --------
    Returns
    --------
    A list containing the sizes of the avalanches and a list containing the durations of the avalanches.
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
        end.append(len(prova))
            
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
        end.append(len(prova))
            
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
    color = cm.Reds_r(np.linspace(0.8,1,10))

    x= np.arange(0,sample.shape[1])
    for i in range(len(times)-1):
  
        if s[i]:
            plt.plot([times[i] for r in range(len(x))],x, color = 'red', linewidth = 0.000001)
            plt.fill_betweenx(x,np.asarray([times[i] for r in range(len(x))]),np.asarray([times[i+1] for r in range(len(x))]), color = color[0])
        #else:
            #plt.plot([times[i] for r in range(len(x))],x, '-', color = 'white', alpha = 0.3)
            #plt.fill_betweenx(x,np.asarray([times[i] for r in range(len(x))]),np.asarray([times[i+1] for r in range(len(x))]), color = 'white')

    
    for j in range(sample.shape[1]):
        idx = np.where(sample[:,j]> 0)[0]
        plt.plot(idx,[j for i in range(len(idx))], 'k|', markersize = 1.5)
        
    plt.xlabel('Time (temporal frame)')
    plt.ylabel('Unit')
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


def scaling(sizes, durations, avinterval,  lim1 = 4 , lim2 = 4,ax = None, tau = "default", errtau = "default", alpha = "default", erralpha = "default", maxxminsizes = "default", maxxmindur = "default", xmaxsizes = "default", xmaxdur = "default"):
    
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
        ax.fill_between(x,(10**inter)*x**(fit + errfit*-3),
                         (10**inter)*x**(fit+ errfit*3),color= greens[3])


        ax.fill_between(x,(10**inter)*x**(pred+errpred*-3),
                         (10**inter)*x**(pred+errpred*+3),color= "pink", alpha = 0.5)
    
    else:
        ax.fill_between(x,(10**inter)*x**(pred+errpred*-3),
                         (10**inter)*x**(pred+errpred*+3),color= "pink", alpha = 0.5)
        
        
        ax.fill_between(x,(10**inter)*x**(fit + errfit*-3),
                         (10**inter)*x**(fit+ errfit*3),color= greens[3])

    
    x_annot_tau = 0.6
    y_annot_tau = 0.15
    ax.annotate(r'$\delta_{fit} =%2.2f \pm %2.2f $'%(round(fit,2), round(errfit,2) ),xy=(x_annot_tau,y_annot_tau),fontsize = 'x-large',xycoords = 'axes fraction')
    
    x_annot_tau = 0.6
    y_annot_tau = 0.25
    
    ax.annotate(r'$\delta_{pred} =%2.2f \pm %2.2f$'%(round(pred,2) ,round(errpred,2) ),xy=(x_annot_tau,y_annot_tau),fontsize = 'x-large',xycoords = 'axes fraction')
    ax.set_xlabel('Durations [# of time bins]')
    ax.set_ylabel('Sizes [# of events]')


    ax.legend(loc = (0,0.6))

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

def picchibinned2(sample1, interv): 
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

    n = np.sum(newvec, axis = 1)
    
    
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

    
def GaussianComparison(sample,thres,coef, n1, n2 = 0):
    """
    Compares the distribution of the amplitudes in one electrode with the best Fit Gaussian
    
    Parameters
    --------
    sample : tri/bidimensional array of the continuous data. Dimensions : temporal dim x spatial dim1 (x spatial dim2)
    Thres : tri/bidimensional array of the standard deviations of the data/medians of the data... Dimensions : temporal dim x spatial dim1 (x spatial dim2)
    coef : multiplicative coefficient for the thresholds
    (n1,n2) : coordinates of the chosen electrode
    
    
    """
    if thres == "stds":
        Thres = np.std(sample,axis = 0,ddof = 1)
    if thres == "median":
         Thres = np.median(sample,axis = 0)
            
    if sample.ndim < 3:
        sample = sample.reshape(sample.shape[0],sample.shape[1],1)
    if Thres.ndim < 2:
        Thres = Thres.reshape(Thres.shape[0],1)
        
    if sample.ndim == 1:
        sample = sample.reshape(sample.shape[0],1,1)
        
    if Thres.ndim ==1:
        Thres = Thres.reshape(1,1)
        
    fig = plt.figure()

    
    plt.title('Electrode in coordinates (%d,%d) [Array dimensions: 4x55] ' %(n1+1,n2+1))

    sig = (sample[:,n1,n2] - np.mean(sample[:,n1,n2]))/(Thres[n1,n2]) # normalizes the signal by the SD

    bins = np.linspace(min(sig),max(sig),300) #attenzione
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
            bins = np.linspace(xlim1,xlim2, 300) #attenzione
            a,b = np.histogram(sig, bins =bins , density = True)
            pdf.append(a)
           
    coord = (b[:-1]+ b[1:])/2 
    """
    plt.plot(coord, pdf[0], '.',color = 'gray', alpha = 0.5, label='Individual channels')
    for s in range(1,len(pdf)): #attenzione
        plt.plot(coord, pdf[s], '.',color = 'gray', alpha = 0.5)
    """
    pdfmean = np.mean(pdf,axis = 0)
    pdfstd = np.std(pdf,axis = 0)/np.sqrt(len(pdf))

    x = coord
    y = pdfmean


    popt,pcov = curve_fit(gaus,x,y,p0=[0,1],bounds = (np.array([0-0.5, -4]),np.array([0+0.5,4])))



    plt.plot(x,gaus(x,*popt),'--', color = 'purple',label='Best Gaussian Fit')


    plt.xlabel('Normalized Amplitudes (SD)')
    plt.ylabel('Probability Density')
    v = np.arange(0,0.6,0.01)

    x1 = np.array([-thres for i in range(len(v))])
    x2 = np.array([thres for i in range(len(v))])
    #plt.plot(x1,v, 'r-', label = r'$\pm$ %s SD' %thres)
    #plt.plot(x2,v, 'r-', r'%s SD' %thres)
    plt.fill_between(x,pdfmean + -3*pdfstd, pdfmean+3*pdfstd, color = 'blue', alpha = 0.3)
    plt.plot(x,pdfmean, '-', color = 'blue', linewidth = 2, label = 'Average')
    
    plt.yscale('log')
    plt.xticks(np.arange(xlim1,xlim2+1,1))
    plt.legend(loc = (1,0.6))

def gaus(x,x0,sigm):
    return 1/(np.sqrt(2*np.pi*sigm**2))*np.e**(-(x-x0)**2/(2*sigm**2))


  
    
