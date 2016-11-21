import matplotlib.pyplot as plt
import numpy as np
from obspy.noise.correlation_functions import phase_xcorr

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
    corr = phase_xcorr(data[0], data[1], maxlag)  

    if plot:
        plt.subplot(211)
        plt.plot(np.arange(len(corr[0])) * 0.05, np.abs(corr[0]))
        plt.subplot(212)
        plt.plot(np.arange(len(corr[1])) * 0.05, np.abs(corr[1]))


    if plot:
        plt.figure()
        plt.plot(corr)


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
