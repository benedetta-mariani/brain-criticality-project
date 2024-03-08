import numpy as np
import matplotlib.pyplot as plt
import powerlaw as pwl
from matplotlib import cm
from statsmodels.regression import linear_model as sm
from scipy.signal import find_peaks
from scipy import signal
from scipy.stats import percentileofscore

def threshold2(sample1,means,stds,thres):
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
            init = np.where(np.diff(prova)>0)[0]
            end = np.where(np.diff(prova)<0)[0]
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


def threshold3(sample1,means,stds,thres):
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
            

            initsign = np.where((changesign)>0)[0]
            endsign = np.where((changesign)<0)[0]

            init = np.where(np.diff(prova)>0)[0]
            end = np.where(np.diff(prova)<0)[0]
            #print(len(initsign) ==len(endsign))
            #print(len(sig),len(initsign),len(init), len(end))
            if len(init) < len(end):
                init = np.insert(init,0,0)

            if len(end) < len(init):
                end = np.append(end,len(sig))

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
            while f < len(end):
                if not len((set(np.arange(init[a], end[f],1)) & set(initsign[:]))):
                    k+=1
                    f+=1

                else:
                    if k> 0:
                        endss.append(end[f-1])
                        initss.append(init[a])
                        a += k
                        k = 0
                        #f+=1
                    else:
                        
                        endss.append(end[f-1])
                        initss.append(init[a])
                        a += 1
                        #f+=1  
            for l in range(len(initss)):
                groups.append((sig[initss[l]:endss[l]]-means[s]))##
                times.append(tempi[initss[l]:endss[l]])
            zeta = []
            for m in range(len(groups)):
                zeta.append(times[m][groups[m].tolist().index(min(groups[m]))])
            sample2[zeta,s] = 1
            
    return sample2.reshape(initshape)
    
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
    return singledur,medie,std

def delta(alpha, salpha, tau, stau):
    """
    returns delta = (alpha - 1)/(tau - 1) and corresponding error estimate
    """
    return (alpha - 1)/(tau -1), np.sqrt((1/(tau -1))**2*salpha**2+ ((1-alpha)/(tau -1)**2)**2*stau**2)


def scaling(sizes, durations, avinterval,  lim1 = 4 , lim2 = 4,ax = None, tau = "default", errtau = "default", alpha = "default", erralpha = "default", maxxminsizes = "default", maxxmindur = "default", xmaxsizes = "default", xmaxdur = "default", xminfit = 'default',xmaxfit = 'default'):
    
    if ax == None:
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
        tau,errtau = exponent(sizes,maxxmin = maxxminsizes,xmax =xmaxsizes,lim = lim1)
        
    if alpha == "default" and erralpha == "default":
        alpha,erralpha = exponent(durations,maxxmin = maxxmindur,xmax =xmaxdur,lim = lim2)

    pred =  delta(alpha, erralpha, tau, errtau )[0]
    errpred = delta(alpha, erralpha,tau, errtau)[1]


    durations = np.array(durations)*avinterval
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

    ax.set_xscale('log')
    ax.set_yscale('log')
   
    x = np.arange(xminfit,xmaxfit,1)

    ax.plot(durations, sizes, '.', color = blues[15],alpha = 0.3)
    ax.plot(np.asarray(durations)[nott],np.asarray(sizes)[nott], '.', color = 'gray', alpha = 0.5)
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
    ax.annotate(r'$\delta_{fit} =%2.2f \pm %2.2f $'%(round(fit,2), round(errfit,2) ),xy=(x_annot_tau,y_annot_tau),xycoords = 'axes fraction',fontsize= 17)
    
    x_annot_tau = 0.6
    y_annot_tau = 0.25
    
    ax.annotate(r'$\delta_{pred} =%2.2f \pm %2.2f$'%(round(pred,2) ,round(errpred,2) ),xy=(x_annot_tau,y_annot_tau),xycoords = 'axes fraction',fontsize= 17)

    ax.legend(loc = (0.01,0.6))
    ax.errorbar(a,b,yerr = c, fmt = 'o', color = blues[29],markersize = 3.5,barsabove = False,capsize = 3, elinewidth = 3,capthick = 1)


def return_param(sample,maxxmin = "default",xmax = "default",lim = 4):

    if maxxmin ==  "default":
        maxxmin = max(sample)
    if xmax ==  "default":
        xmax = max(sample)
  
    ypred = pwl.Fit(sample, xmin = (1,maxxmin+1),xmax = xmax, parameter_range = {'alpha' : [1.,lim]}, discrete = True)
    return ypred.power_law.xmin, ypred.power_law.alpha
    

def exponent(sample,maxxmin = "default",xmax = "default",lim = "default"):
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
    xmin here is the specific value that has to be used as xmin
    """
    ypred = pwl.Fit(sample,xmin = (xmin,xmin + 1),xmax = xmax,parameter_range = {'alpha' : [1,lim]},discrete = True)
    return ypred

def goodness(sample, xmin, xmax, Ds,lim):
    """
    given N randomly generated Ds (syntheticly generated Kolmogorov Smirov distances) computes the goodness of the power law fit
    """
    predic = ypred(sample, xmin, xmax,lim)
    Dreal = predic.power_law.D
    score = percentileofscore(Ds,Dreal)
    pval = 1-score/100.
    return pval

def modelcompare(sample,xmin,xmax,model1,model2):
    ypred = pwl.Fit(sample,xmin =(1,xmin+1),xmax = xmax,discrete = True)
    return ypred.distribution_compare(model1, model2, normalized_ratio = True)
def returnbin(n,interv):
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


def RasterPlot(sample, xtime, av, ax = None, av_color = 'red',
               nsize = 1.5, alpha = 0.3):
    """
    Parameters
    --------
    sample : Array of discretized data. Dimensions : temporal dim x spatial dim1 (x spatial dim2)
    av : width of temporal bin 
    Returns
    --------
    Plots the Raster Plots and the detected avalanches (an avalanche is preceded and followed by white bins)
    """
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    s = returnbin(events(sample),av)
    if s[0] == 1:
        start_av = True
    else:
        start_av = False
    s = xtime[(np.where(np.diff(s) != 0)[0]+1)*av]
    if start_av:
        s = np.concatenate([[0], s])
    if s.size % 2 != 0:
        s = np.concatenate([s, [s[-1] + av]])

    for sval in s.reshape(-1, 2):
        ax.axvspan(sval[0], sval[1], 0.01, 0.99, color = av_color, alpha = alpha, zorder = -10)
    for j in range(sample.shape[1]):
        idx = np.where(sample[:,j]> 0)[0]
        ax.plot(xtime[idx], np.ones(len(idx))*j, '|', color = 'black', markersize = nsize)

def KuramotoIndex(sample, nunits = "default", time = "default"):
    """
    returns Kuramoto parameter from continuous data
    sample.shape = time x nunits
    """
    if sample.ndim > 2:
        sample = sample.reshape(-1, sample.shape[1]*sample.shape[2])
    if nunits == "default":
        nunits = sample.shape[1]
    if time == "default":
        time = sample.shape[0]
    
    sig = signal.hilbert(sample, axis = 0)
    phases = np.angle(sig)
    sincr = np.abs(np.mean(np.exp(phases*1j),1)) #Kuramoto indexes at each time
    return np.mean(sincr),sincr

def calc_cv(SpikeMat):
    '''Calculate the Inter-Spike Intervals (ISI) and the
    Coefficient of Variation (CV) from given spike timings.
    SpikeMat = numpy.array, shape = (time,nunits)
    
    cv = 1: poisson process-like irregularity. std of ISI is of the same order of magnitude of the mean. irregular regime
    cv = 0: completely regular regime. std of ISIs is 0.
    cv < 1: more regular than a poisson process.
    cv > 1: high variability in ISIs, more irregular than a poisson process.

    '''

    if SpikeMat.shape[1] > SpikeMat.shape[0]:
        raise Exception('Error, the array must be transposed (first dimension should be time)')

    cvs =[]
    isis = []
    for nrn in range(SpikeMat.shape[1]):
        spikes = np.where(SpikeMat[:,nrn]>0)[0]
        times = np.diff(spikes)
        isis += times.tolist()
        CV = np.std(times)/np.mean(times)
        cvs.append(CV)
    return cvs,isis

def get_phase(SpikeMat, time_array):
    ''' 
    Get the phase values per neuron based on the spike timings
    SpikeMat = numpy.array, shape = (time, # of units)
    time_array = array of time, e.g. np.arange(0,T,dt) or np.arange(0,N,1)
    '''
    if SpikeMat.shape[1] > SpikeMat.shape[0]:
        raise Exception('Error, the array must be transposed (first dimension should be time)')

    nunits = min(SpikeMat.shape) # number of neurons
    phase = np.empty((nunits, len(time_array)))

    phase[:] = np.nan
    for nrn in np.arange(0,nunits):
        tspike = time_array[np.array(SpikeMat[:,nrn],dtype =bool)]
        ispike = np.where(SpikeMat[:,nrn])
        
        for ii in np.arange(1,len(ispike[0])):
            idx1 = ispike[0][ii-1]
            idx2 = ispike[0][ii]
            tv_tmp = (time_array[idx1:idx2]-time_array[idx1])/(time_array[idx2]-time_array[idx1]) # normalize in [0,1]
            phase[nrn,idx1:idx2] = 2*np.pi*tv_tmp # normalize in [0,2π]
        
    return phase
    
def Kuramoto_param(phases):
    """
    phases shape = (nunits,time)
    """
    if sample1.shape[0] > sample1.shape[1]:
        raise Exception('Error, the array must be transposed (time should be the second dimension)')
    r = np.abs(np.nanmean(np.exp(1j*phases),axis=0))
    return r
