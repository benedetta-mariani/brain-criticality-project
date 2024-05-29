""" 
more complete set of functions wrt avalanches.py to study avalanches statistics 
(e. g. comprehend study of avalanches profiles and calculation of inter avlanches times)
"""

import numpy as np
import matplotlib.pyplot as plt
import powerlaw as pwl
from matplotlib import cm
from statsmodels.regression import linear_model as sm
from scipy.signal import find_peaks
from scipy import signal
from scipy.stats import percentileofscore

def compute_avalanches(data,interv,dt =1):
    final_tt = data.T
    ev = np.sum(final_tt,0)

    if len(ev)%interv > 0:
        add = (int(len(ev)/interv) + 1)* interv - len(ev)
        ev = ev.tolist()
        for i in range(add):
            ev = ev + [0]

    ev = np.asarray(ev).reshape(int(len(ev)/interv), interv)
    global_signal =np.sum(ev,axis = 1)
    new= np.array(global_signal, dtype = bool)
    ##ev.shape, new.shape, mean_interspike_time
    #print(np.array(new, dtype= int)[:15])
    final_t = np.array(new, dtype = float)
    
    av_indice_start = np.where((final_t[1:] - final_t[:-1])>0)[0]  + 1# These are the indices where an avalanche begins
    av_indice_end = np.where((final_t[1:] - final_t[:-1])< 0)[0] + 1  # hese are the indices where an avalanche ends
    
    #We modify the content of both av_indice_start and av_indice_end depending on the 
    #initial and final values of final_t
    
    #Here both av_indice_start and av_indice_end register correctly all avalanche begginings and ends
    if final_t[0]==0 and final_t[-1]==0: 
        av_indice_start = av_indice_start
        av_indice_end = av_indice_end  # No changes on both arrays since they are correct at possitions and lenghts 
    
    #Here av_indice_end does not register the position of the end of the last avalanche  
    elif final_t[0]==0 and final_t[-1]==1:
        av_indice_end = np.append(av_indice_end,len(final_t)-1) 
    
    #Here av_indice_start does not register the position of the begining of the first avalanche
    elif final_t[0]==1 and final_t[-1]==0:
        av_indice_start = np.insert(av_indice_start,0,0)
    
    #Here av_indice_start does not register the position of the begining of the first avalanche
    # and av_indice_end does not register the position of the end of the last avalanche
    
    elif final_t[0]==1 and final_t[-1]==1:
        av_indice_start = np.insert(av_indice_start,0,0)
        av_indice_end = np.append(av_indice_end,len(final_t)-1)
               
        
   
    t=np.arange(0,len(final_t)*interv/25000,0.00004)
    print(t[-1])
    avalanche_sizes = []
    avalanche_dur = []

    for s in range(len(av_indice_start)):
        if len(av_indice_start) != len(av_indice_end):
            print('Error, they must be of the same length')
            break
        avalanche_sizes.append(np.sum(global_signal[av_indice_start[s]:av_indice_end[s]]))
        avalanche_dur.append((t[av_indice_end[s]]- t[av_indice_start[s]])/dt)##
    return np.array(avalanche_sizes),np.array(avalanche_dur)

def intertimes(data,interv,dt =1):
    final_tt = data.T
    ev = np.array(np.sum(final_tt,0), dtype=bool)
    if len(ev)%interv > 0:

        add = (int(len(ev)/interv) + 1)* interv - len(ev)
        ev = ev.tolist()
        for i in range(add):
            ev = ev + [0]

    ev = np.asarray(ev).reshape(int(len(ev)/interv), interv)
    new = np.array(np.sum(ev,axis = 1),dtype = bool)
    ##ev.shape, new.shape, mean_interspike_time
    #print(np.array(new, dtype= int)[:15])
    final_t = np.array(new, dtype = float)
    
    #The key part
    av_indice_start = np.where((final_t[1:] - final_t[:-1]) <0)[0]  + 1# These are the indices where an avalanche begins
    av_indice_end = np.where((final_t[1:] - final_t[:-1])> 0)[0] + 1  # hese are the indices where an avalanche ends
    
    #We modify the content of both av_indice_start and av_indice_end depending on the 
    #initial and final values of final_t
    
    #Here both av_indice_start and av_indice_end register correctly all avalanche begginings and ends
    if final_t[0]==1 and final_t[-1]==1: 
        av_indice_start = av_indice_start
        av_indice_end = av_indice_end  # No changes on both arrays since they are correct at possitions and lenghts 
    
    #Here av_indice_end does not register the position of the end of the last avalanche  
    elif final_t[0]==1 and final_t[-1]==0:
        av_indice_end = np.append(av_indice_end,len(final_t)-1) 
    
    #Here av_indice_start does not register the position of the begining of the first avalanche
    elif final_t[0]==0 and final_t[-1]==1:
        av_indice_start = np.insert(av_indice_start,0,0)
    
    #Here av_indice_start does not register the position of the begining of the first avalanche
    # and av_indice_end does not register the position of the end of the last avalanche
    
    elif final_t[0]==0 and final_t[-1]==0:
        av_indice_start = np.insert(av_indice_start,0,0)
        av_indice_end = np.append(av_indice_end,len(final_t)-1)
               
        
   
    t=np.arange(0,len(final_t)*interv/25000,0.00004)
    #print(t[-1])
    #avalanche_sizes = []
    avalanche_dur = []

    for s in range(len(av_indice_start)):
        if len(av_indice_start) != len(av_indice_end):
            print('Error, they must be of the same length')
            break
        #avalanche_sizes.append(np.sum(global_signal[av_indice_start[s]:av_indice_end[s]]))
        avalanche_dur.append((t[av_indice_end[s]]- t[av_indice_start[s]])/dt)##??? giusto così
        
    return np.array(avalanche_dur)

def thresholdnuova(sample1,means,stds,thres):
    """ 
    Detects as events the points of maximum excursion over a threshold, considering either positive and negative excursions or only negative. if "option1" is selected, the one largest maximum between two crossings of the mean assigns the final event time.
    For a faster thresholding use the function below findpeaks.
    
    Parameters
    --------
    sample1 : tri or bidimensional array of recorded voltages (or even single time series from one electrode) (shape = temporal dim x spatial dim1 (x  spatial dim2) ) 
    means : array of the means of the signals (shape = spatial dim1 (x spatial dim2) ) 
    stds : array of the thresholds for each channel (standard deviations/medians...) (shape = spatial dim1 (x spatial dim2)
    thres : multiplicative coefficent for the thresholds
    choose : if "posneg" both positive and negative deflections are considered, if "neg" only negative
    opz : if "option1" the one largest maximum between two crossings of the mean assigns the final event time, if "option2" an event is simply the point of maximum excursion over a threshold
    
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
        
    if sample1.shape[1] > sample1.shape[0]:
         raise Exception('Error, the array must be transposed (first dimension should be time)')

    sample2 = np.zeros(sample1.shape, dtype = int)

    for s in range(sample1.shape[1]):
        if stds[s]>0:
            sig = sample1[:,s].reshape(-1)
            stan = stds[s]
            tempi = np.arange(0,len(sig),1)
            prova =np.array((sig - means[s]) <= -thres*stan, dtype = float)
            init = np.where(np.diff(prova)>0)[0] + 1
            end = np.where(np.diff(prova)<0)[0] + 1
            if len(init) < len(end):
                init = np.insert(init,0,0)

            if len(end) < len(init):
                end = np.append(end,len(sig))

            #print(len(init), len(end))
            groups = []
            times = []
            for l in range(len(init)):
                groups.append(np.abs(sig[init[l]:end[l]]-means[s]))##
                times.append(tempi[init[l]:end[l]])
            #print(groups)
            zeta = []
            for m in range(len(groups)):
                zeta.append(times[m][groups[m].tolist().index(max(groups[m]))])
            sample2[zeta,s] = 1
            
    return sample2.reshape(initshape)

def thresholdnuova2(sample1,means,stds,thres):
    """ 
    Detects as events the points of maximum excursion over a threshold, considering either positive and negative excursions or only negative. if "option1" is selected, the one largest maximum between two crossings of the mean assigns the final event time.
    For a faster thresholding use the function below findpeaks.
    
    Parameters
    --------
    sample1 : tri or bidimensional array of recorded voltages (or even single time series from one electrode) (shape = temporal dim x spatial dim1 (x  spatial dim2) ) 
    means : array of the means of the signals (shape = spatial dim1 (x spatial dim2) ) 
    stds : array of the thresholds for each channel (standard deviations/medians...) (shape = spatial dim1 (x spatial dim2)
    thres : multiplicative coefficent for the thresholds
    choose : if "posneg" both positive and negative deflections are considered, if "neg" only negative
    opz : if "option1" the one largest maximum between two crossings of the mean assigns the final event time, if "option2" an event is simply the point of maximum excursion over a threshold
    
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
        
    if sample1.shape[1] > sample1.shape[0]:
         raise Exception('Error, the array must be transposed (first dimension should be time)')

    sample2 = np.zeros(sample1.shape, dtype = int)

    for s in range(sample1.shape[1]):
        if stds[s]>0:
            sig = sample1[:,s].reshape(-1)
            stan = stds[s]
            tempi = np.arange(0,len(sig),1)
            prova =np.array((sig - means[s]) <= -thres*stan, dtype = float)
            
            changesign = np.diff(np.sign(sig[:]-means[s]))
            
            
            #changesign1 = np.hstack((changesign1,0.))
            

            initsign = np.where((changesign)<0)[0] + 1
            endsign = np.where((changesign)>0)[0] + 1

            init = np.where(np.diff(prova)>0)[0] + 1
            end = np.where(np.diff(prova)<0)[0] + 1
            #print(len(initsign) ==len(endsign))
            #print(len(sig),len(initsign),len(init), len(end))
            if len(init) < len(end):
                init = np.insert(init,0,0)

            if len(end) < len(init):
                end = np.append(end,len(sig)) ###???

            #print(init, end)
            groups = []
            times = []
            
            #intss = []
            #endss = []
            
                        
            initss = []
            endss = []

            a = 0
            f = 0
            k = 0
            endss = []
            initss = []
            #g = 0 
            o = 0
            while f < len(end):
                if not len((set(np.arange(init[a], end[f],1)) & set(initsign[:]))):
                    if (init[a]> end[f]):print('STRANGe')
                    
                    k+=1
                    f+=1
                   # o+=1

                else:
                    if k> 0:
                        endss.append(end[f-1])
                        initss.append(init[a])
                        a += (k)
                        k = 0
                
            if not len((set(np.arange(init[a], end[-1],1)) & set(initsign[:]))):
                        endss.append(end[-1])
                        initss.append(init[a])
     
            for l in range(len(initss)):
                groups.append((sig[initss[l]:endss[l]]-means[s]))##
                times.append(tempi[initss[l]:endss[l]])
            zeta = []
            for m in range(len(groups)):
                zeta.append(times[m][groups[m].tolist().index(min(groups[m]))])
            sample2[zeta,s] = 1
            
    return sample2.reshape(initshape)


def findpeaks(sig, thres, choose, dist = 200):
    """
    Finds peaks in a time series.
    Set eventual other constraints in the scipy.signal function find_peaks
    """
    sig2 = np.zeros(sig.shape, dtype = int)
    if choose == "neg":
        p1 = find_peaks(-(sig - np.mean(sig)), height = thres, distance = dist)[0]
        if np.any(p1):
            sig2[p1] = 1
            
    if choose == "posneg":
        p1 = find_peaks(-(sig - np.mean(sig)), height = thres, distance = dist)[0]
        p2 = find_peaks((sig - np.mean(sig)), height = thres, distance = dist)[0]
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
    if sample.shape[1] > sample.shape[0]:
         raise Exception('Error, the array must be transposed (first dimension should be time')
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
    if n.ndim > 1:
        raise Exception('Error! array must be one-dimensional')
    else:
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
    prova = np.array(new>0, dtype = int)
 
    start = np.where(np.diff(prova)>0)[0] + 1 
    end = np.where(np.diff(prova)<0)[0] +1
    # print(len(new), len(start),len(np.diff(prova)), len(prova))
    # print(start[:10])
    
    if prova[0] == 1:
        start = np.insert(start,0,0)
    if prova[-1] == 1:
        end = np.append(end,len(new)) ## check
    
    #print(end[-10:])
    avalanches = []
    for s in range(len(start)):
        avalanches.append(new[start[s]:end[s]])
        #print(new[start[s]:end[s]])
    sizes = []
    durations = []
    for l in range(len(avalanches)):
        sizes.append(np.sum(avalanches[l]))
        durations.append(len(avalanches[l]))
    
    return sizes, durations

    

    
def main(sample,av = "default"):
    """
    Parameters
    --------
    sample : discretized array of events
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


def intertempi(n):
    if n.ndim > 1:
        raise Exception('Error! array must be one-dimensional')
    idx = np.where(n > 0)[0]
    intertempi = idx[1:] - idx[:-1] 
    return intertempi

def sgivent(sizes, durations):
    """
    returns unique values of durations and corresponding average sizes of avalanches
    """

    tot = [durations,sizes]

    singledur = np.sort(np.unique(durations))
    sdit = [[] for i in range(len(singledur))]
    for i in range(len(singledur)):
        for z in range(len(durations)):
            if tot[0][z] == singledur[i]:
                sdit[i].append(tot[1][z])

    medie = []
    std = []
    for i in range(len(sdit)):
        medie.append(np.mean(sdit[i]))
        std.append(np.std(sdit[i])/np.sqrt(len(sdit[i])))
    return singledur,np.array(medie),np.array(std)

def delta(alpha, salpha, tau, stau):
    """
    returns delta = (alpha - 1)/(tau - 1) and corresponding error estimate
    """
    return (alpha - 1)/(tau -1), np.sqrt((1/(tau -1))**2*salpha**2+ ((1-alpha)/(tau -1)**2)**2*stau**2)


def scaling(sizes, durations,  lim1 = 4 , lim2 = 4,ax = None, tau = "default", errtau = "default", alpha = "default", erralpha = "default", maxxminsizes = 100, maxxmindur = 0.05, xmaxsizes = "default", xmaxdur = "default", xminfit = 'default',xmaxfit = 'default', plotto = False):
    
    if ax == None:
        if plotto:
            fig = plt.figure(figsize = (6,4))
            ax = fig.add_subplot(1,1,1)
        
    if maxxminsizes ==  "default":
        maxxminsizes = max(sizes)
    if xmaxsizes ==  "default":
        xmaxsizes = max(sizes)

    
    if maxxmindur ==  "default":
        maxxmindur = max(durations)
    if xmaxdur ==  "default":
        xmaxdur = max(durations)
    
    
    if tau == "default" and errtau == "default":
        tau,errtau = exponent(sizes,maxxmin = maxxminsizes,xmax =xmaxsizes,lim = lim1,)
        
    if alpha == "default" and erralpha == "default":
        alpha,erralpha = exponent(durations,maxxmin = maxxmindur,xmax =xmaxdur,lim = lim2)

    pred =  delta(alpha, erralpha, tau, errtau )[0]
    errpred = delta(alpha, erralpha,tau, errtau)[1]


    durations = np.asarray(durations)
    if xminfit== 'default': xminfit = min(durations)
    if xmaxfit== 'default': xmaxfit = max(durations)
    
    
    prova = np.array([np.asarray(sizes), np.asarray(durations)])
    prova = prova.transpose()
    prova2 = [0 for i in range(len(prova))]
    
    for r in range(len(prova)):
        if prova[r][1] < xminfit  or prova[r][1] > xmaxfit :
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
        if prova[r][1] < xminfit or prova[r][1] > xmaxfit :

            nott[r] = True
        else:

            nott[r] = False

    print('Prediction from crackling noise relation: delta = ',pred, '+-', errpred)
    print('Fit from of average size given duration points: delta = ',fit, '+-', errfit)

   
    x = durations
    #plotto = False
    #np.arange(xminfit,xmaxfit,0.002)
    if plotto:

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.plot(durations, sizes, '.', color = blues[15],alpha = 0.3)
        ax.plot(np.asarray(durations)[nott],np.asarray(sizes)[nott], '.', color = 'gray', alpha = 0.5)
        ax.plot(x, (10**inter)*x**pred, 'r', label = 'Prediction', lw = 2)
        ax.plot(x, (10**inter)*x**fit, color = greens[12], label = 'Fit', lw = 2)
        #if errfit > errpred:
        #    ax.fill_between(x,(10**inter)*x**(fit + errfit*-3),
        #                     (10**inter)*x**(fit+ errfit*3),color= greens[3])


    #        ax.fill_between(x,(10**inter)*x**(pred+errpred*-3),
                             #(10**inter)*x**(pred+errpred*+3),color= "pink", alpha = 0.5)

       # else:
       #     ax.fill_between(x,(10**inter)*x**(pred+errpred*-3),
                             #(10**inter)*x**(pred+errpred*+3),color= "pink", alpha = 0.5)


          #  ax.fill_between(x,(10**inter)*x**(fit + errfit*-3),
                             #(10**inter)*x**(fit+ errfit*3),color= greens[3])

        x_annot_tau = 0.6
        y_annot_tau = 0.15
        ax.annotate(r'$\delta_{fit} =%2.2f \pm %2.2f $'%(round(fit,2), round(errfit,2) ),xy=(x_annot_tau,y_annot_tau),xycoords = 'axes fraction',fontsize= 17)

        x_annot_tau = 0.6
        y_annot_tau = 0.25

        ax.annotate(r'$\delta_{pred} =%2.2f \pm %2.2f$'%(round(pred,2) ,round(errpred,2) ),xy=(x_annot_tau,y_annot_tau),xycoords = 'axes fraction',fontsize= 17)
        ax.set_xlabel('Avalanche durations')
        ax.set_ylabel('Avalanche sizes')
        
        ax.legend(loc = (0.01,0.6))
        ax.errorbar(a,b,yerr = c, fmt = 'o', color = blues[29],markersize = 3.5,barsabove = False,capsize = 3, elinewidth = 3,capthick = 1)
    return np.asarray([pred, errpred,fit,errfit,a,b,c])


def RasterPlot(sample, av, ticks,ax, color, alpha):
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
    #color = cm.Reds_r(np.linspace(0.8,1,10))

    x= np.arange(0,sample.shape[1])
    for i in range(len(times)-1):
  
        if s[i]:
            plt.plot([times[i] for r in range(len(x))],x, color = 'red', linewidth = 0.000001)
            plt.fill_betweenx(x,np.asarray([times[i] for r in range(len(x))]),np.asarray([times[i+1] for r in range(len(x))]), color = color, alpha = alpha)
        #else:
            #plt.plot([times[i] for r in range(len(x))],x, '-', color = 'white', alpha = 0.3)
            #plt.fill_betweenx(x,np.asarray([times[i] for r in range(len(x))]),np.asarray([times[i+1] for r in range(len(x))]), color = 'white')

    
    for j in range(sample.shape[1]):
        idx = np.where(sample[:,j]> 0)[0]
        plt.plot(idx,[j for i in range(len(idx))], '|',color = black, markersize = 1.5)
        
    #plt.xlabel('Temporal frame')
    plt.ylabel(r'Unit')
    ax.spines['bottom'].set_color('None')
    plt.xticks(ticks)
    plt.yticks()
    
    


def avalanche_finder(S_shape_, coef):
    where_spikes = np.where(S_shape_ != 0)
    interspike_time = (where_spikes - np.roll(where_spikes,1))
    interspike_time = np.delete(interspike_time,0) # remove the first element
    mean_interspike_time = np.sum(interspike_time)/len(interspike_time)
    mean_interspike_time = int(int(round(mean_interspike_time))*coef)
    n = len(S_shape_)
    # Now I modify the S_shape_ in order to applay the ney binning
    rest = int(n%mean_interspike_time)
    if rest == 0:
        n_reduced = int(n/mean_interspike_time)
    else:
        n_reduced = int((n-rest)/mean_interspike_time)
    S_shape = np.zeros(n_reduced)
    delta = mean_interspike_time
    for i in range(0, n-rest, mean_interspike_time):
        j = int(i/mean_interspike_time)
        S_shape[j] = S_shape_[i:i+delta].sum()
    t_in = [] # timestep in which avalanche starts
    t_fin = [] #                            ends
    avalanche = False
    for timestep in range(n_reduced):
        if S_shape[timestep] != 0:
            if avalanche == False:
                t_in.append(timestep) # lists faster than array in this
                avalanche = True
        if S_shape[timestep] == 0 and avalanche == True:
            t_fin.append( timestep ) 
            avalanche = False
    if len(t_in) != len(t_fin):
        t_in.pop(-1)
        # convert them into arrays for efficiency purposes
    t_in = np.array(t_in)
    t_fin = np.array(t_fin)
    sizes = []
    for i in range(len(t_in)):
        if S_shape[t_in[i]:t_fin[i]].sum() == 0: print("ERROR: zeros inside sizes!") # control
        sizes.append( S_shape[t_in[i]:t_fin[i]].sum() )
    durations = np.array(t_fin - t_in)  # number of bins
    duration_uni = np.unique(durations).astype(int) # uniques durations of avalanches
    shape_mean = [] # mean temporal shapes
    frequencies = []
    dur = duration_uni[:]
    shape_mean = [[] for g in range(len(dur))]
    e = 0
    for d in dur:
        shape_temp = np.zeros(d)
        counter = 0
        for i in range(len(t_in)):
            if t_fin[i] - t_in[i] == d:
                shape_temp = S_shape[t_in[i]:t_fin[i]]
                counter += 1
                shape_mean[e].append(shape_temp)
                
        frequencies.append(counter)
        e += 1
    shape_mean = np.array(shape_mean)
    return sizes, durations, S_shape,shape_mean, frequencies

def exponent(sample,maxxmin = "default",xmax = "default",lim = "default"):
    if maxxmin ==  "default":
        maxxmin = max(sample)
    if xmax ==  "default":
        xmax = max(sample)
    if lim ==  "default":
        lim = 4
        
    ypred = pwl.Fit(sample,xmin = (1,maxxmin + 1),xmax = xmax, parameter_range = {'alpha': [1,lim]}, discrete = True)
    return ypred.power_law.alpha, ypred.power_law.sigma
