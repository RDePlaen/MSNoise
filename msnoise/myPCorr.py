import matplotlib.pyplot as plt
import numpy as np
from obspy.core import Trace, Stream, read
import subprocess
import os



def myPCorr(data, maxlag, plot=False):
    """This function takes ndimensional *data* array, computes the cross-correlation in the frequency domain
    and returns the cross-correlation function between [-*maxlag*:*maxlag*].

    :type data: :class:`numpy.ndarray`
    :param data: This array contains the fft of each timeseries to be cross-correlated.
    :type maxlag: int
    :param maxlag: This number defines the number of samples (N=2*maxlag + 1) of the CCF that will be returned.

    :rtype: :class:`numpy.ndarray`
    :returns: The cross-correlation function between [-maxlag:maxlag]
    """

    stats = {'sampling_rate': 20.}

    st1 = Stream(traces = [Trace(data = data[0], header=stats)])
    st1.write('trace1.sac', format = 'SAC')
 
    st2 = Stream(traces = [Trace(data = data[1], header=stats)])
    st2.write('trace2.sac', format = 'SAC')
   
    maxlag = np.round(maxlag)
    subprocess.call(['Corr_stack_v03.5/pcc5g trace1.sac trace1.sac isac osac pcc tl1=-120. tl2=120. osac'], shell = True)
    pcc_st = read('pcc.sac')
    corr = pcc_st[0].data    
    

    if plot:
        plt.subplot(211)
        plt.plot(np.arange(len(corr[0])) * 0.05, np.abs(corr[0]))
        plt.subplot(212)
        plt.plot(np.arange(len(corr[1])) * 0.05, np.abs(corr[1]))

    if plot:
        plt.figure()
        plt.plot(corr)

    if os.path.isfile('trace1.sac'):
        os.remove('trace1.sac')
    if os.path.isfile('trace2.sac'):
        os.remove('trace2.sac')
    if os.path.isfile('pcc.sac'):
        os.remove('pcc.sac')     
    del data
    return corr


if __name__ == "__main__":
    import time

    data = np.random.random((2, 86400 * 20))
    print(data.shape)
    t = time.time()
    corr = myPCorr(data, maxlag=25, plot=False)
    print("Time:", time.time() - t)
    print(np.mean(corr))

    plt.figure()
    plt.plot(corr)
    plt.axhline(1)
    plt.show()
