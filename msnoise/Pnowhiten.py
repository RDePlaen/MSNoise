import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
from obspy.core import Trace
from scipy import signal

from api import nextpow2

def nowhiten(data, Nfft, delta, freqmin, freqmax, plot=False):
    """This function takes 1-dimensional *data* timeseries array
    and returns the fft between *freqmin* and *freqmax*.

    :type data: :class:`numpy.ndarray`
    :param data: Contains the 1D time series to whiten
    :type Nfft: int
    :para Nfft: The number of points to compute the FFT
    :type delta: float
    :param delta: The sampling frequency of the `data`
    :type freqmin: float
    :param freqmin: The lower frequency bound
    :type freqmax: float
    :param freqmax: The upper frequency bound
    :type plot: bool
    :param plot: Whether to show a raw plot of the action (default: False)

    :rtype: :class:`numpy.ndarray`
    :returns: The FFT of the input trace between the frequency bounds
"""

    if plot:
        plt.subplot(211)
        plt.plot(np.arange(len(data)) * delta, data)
        plt.xlim(0, len(data) * delta)
        plt.title('Input trace')

    stats = {'sampling_rate': 20.}
    trdata=Trace(data, header=stats)
    # Napod = 100
    # Nfft = int(Nfft)
    # freqVec = scipy.fftpack.fftfreq(Nfft,d=delta)[:Nfft/2]
    #
    # J = np.where((freqVec >= freqmin) & (freqVec <= freqmax))[0]
    # low = J[0] - Napod
    # if low <= 0:
    #     low = 1
    #
    # porte1 = J[0]
    # porte2 = J[-1]
    # high = J[-1] + Napod
    # if high > Nfft / 2:
    #     high = Nfft // 2
    trdata.filter("bandpass",freqmin=freqmin,freqmax=freqmax,corners=4,zerophase=True)
#Raph    FFTRawSign = scipy.fftpack.fft(trdata, Nfft)
    
    if plot:
        plt.subplot(212)
        axis = np.arange(len(FFTRawSign))
        plt.plot(axis[1:], np.abs(FFTRawSign[1:]))
        plt.xlim(0, max(axis))
        plt.title('FFTRawSign')

    # # Left tapering:
    # FFTRawSign[0:low] *= 0
    # FFTRawSign[low:porte1] = np.cos(np.linspace(np.pi / 2., np.pi, porte1 - low)) ** 2 * FFTRawSign[low:porte1]
    # # Pass band:
    # FFTRawSign[porte1:porte2] = FFTRawSign[porte1:porte2]
    # # Right tapering:
    # FFTRawSign[porte2:high] = np.cos(np.linspace(0., np.pi / 2., high - porte2)) ** 2 * FFTRawSign[porte2:high]
    # FFTRawSign[high:Nfft+1] *= 0
    
    # Hermitian symmetry (because the input is real)
    #FFTRawSign[-Nfft/2+1:] = FFTRawSign[1:Nfft/2].conjugate()[::-1]

    # if plot:
    #     plt.subplot(413)
    #     axis = np.arange(len(FFTRawSign))
    #     plt.axvline(low, c='g')
    #     plt.axvline(porte1, c='g')
    #     plt.axvline(porte2, c='r')
    #     plt.axvline(high, c='r')
    #
    #     plt.axvline(Nfft - high, c='r')
    #     plt.axvline(Nfft - porte2, c='r')
    #     plt.axvline(Nfft - porte1, c='g')
    #     plt.axvline(Nfft - low, c='g')
    #
    #     plt.plot(axis, np.abs(FFTRawSign))
    #     plt.xlim(0, max(axis))
    #
    #     wdata = np.real(scipy.fftpack.ifft(FFTRawSign))
    #     plt.subplot(414)
    #     plt.plot(np.arange(len(wdata)) * delta, wdata)
    #     plt.xlim(0, len(wdata) * delta)
    #     plt.show()
    
#Raph    return FFTRawSign
    return trdata

    
if __name__ == '__main__':
    import time
    N = 2048
    np.random.seed(1234)
    a = np.random.random(N)
    a = np.sin(a) + np.sin(a / 4.) + np.sin(a / 16.)
    a -= a.mean()
    t = time.clock()
    for i in range(1000):
        nowhiten(a.copy(), N, 0.05, 1.0, 5.9, plot=False)
    print "1000 loops:", (time.clock()-t) * 1000, "ms"
    nowhiten(a.copy(), N, 0.05, 1.0, 5.9, plot=True)
