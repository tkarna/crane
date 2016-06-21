
"""
Computes the spectrum of given signal.

Tuomas Karna 2012-10-04
"""
from crane import plt
import numpy as np
import matplotlib.mlab
from scipy.signal import lombscargle


def computeSpectrum(x, y, freq):
    # compute the spectrum
    pgram = lombscargle(x, y, 2 * np.pi * freq)
    return np.sqrt(4 * pgram / len(x))


def test():
    # First define some input parameters for the signal:

    A = 2.  # amplitude
    w = 1. / 10.0  # sampling frequency (Hz)
    phi = 0.5 * np.pi  # phase
    nin = 1500
    nout = 1000
    frac_points = 0.9  # Fraction of points to discard

    # Randomly select a fraction of an array with timesteps:

    r = np.random.rand(nin)
    x = np.linspace(0.01, 100, nin)
    x_lossy = x[r >= frac_points]

    # Plot a sine wave for the selected times:

    y = A * np.sin(2 * np.pi * w * x + phi)
    y_lossy = A * np.sin(2 * np.pi * w * x_lossy + phi)

    # Define the array of frequencies for which to compute the periodogram:

    f = np.linspace(0.01, 1, nout)

    # Calculate Lomb-Scargle periodogram:

    pg = computeSpectrum(x, y, f)
    pg_lossy = computeSpectrum(x_lossy, y_lossy, f)

    # Now make a plot of the input data:
    pg_fft, f_fft = matplotlib.mlab.psd(
        y, NFFT=16 * 256, Fs=3. / 2 * 1. / w)  # only for evenly sampled data
    plt.subplot(2, 1, 1)
    plt.plot(x, y, 'r.-')
    plt.hold(True)
    plt.plot(x_lossy, y_lossy, 'bo')
    plt.xlabel('time (s)')

    # Then plot the normalized periodogram:

    ax2 = plt.subplot(2, 1, 2)
    #plt.plot(f, sqrt(4*(pg/len(x))))
    #plt.plot(f, pg/len(x))
    ax2.plot(1. / f, pg, 'r')
    plt.hold(True)
    ax2.plot(1. / f, pg_lossy, 'b')
    #ax2.plot(1./f_fft, sqrt(pg_fft),'g')
    ax2.set_xlim([0, 20])
    ax2.legend(['100% data', str(int(round((1 - frac_points) * 100))) + '% data'])
    ax2.set_xlabel('period (s)')
    ax2.set_ylabel('amplitude')
    plt.show()

if __name__ == '__main__':
    test()