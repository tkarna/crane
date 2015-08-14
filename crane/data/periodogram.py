
"""
Simple implementation of a periodogram, following scipy revision 7141
http://projects.scipy.org/scipy/browser/trunk/scipy/signal/spectral.pyx?rev=7141

Tuomas Karna 2012-10-04
"""

from numpy import *

def lombscargle( x, y, freqs) :
  """Computes the Lomb-Scargle periodogram.

  The Lomb-Scargle periodogram was developed by Lomb [1]_ and further
  extended by Scargle [2]_ to find, and test the significance of weak
  periodic signals with uneven temporal sampling.

  The computed periodogram is unnormalized, it takes the value
  ``(A**2) * N/4`` for a harmonic signal with amplitude A for sufficiently
  large N.

  Parameters
  ----------
  x : array_like
      Sample times.
  y : array_like
      Measurement values.
  freqs : array_like
      Angular frequencies for output periodogram.

  Returns
  -------
  pgram : array_like
      Lomb-Scargle periodogram.
  """

  # Check input sizes
  if x.shape[0] != y.shape[0]:
    raise ValueError("Input arrays do not have the same size.")

  # Create empty array for output periodogram
  pgram = zeros(freqs.shape[0], dtype=float64)

  # Local variables
  #cdef Py_ssize_t i, j
  #cdef double c, s, xc, xs, cc, ss, cs
  #cdef double tau, c_tau, s_tau, c_tau2, s_tau2, cs_tau

  for i in range(freqs.shape[0]):

      xc = 0.
      xs = 0.
      cc = 0.
      ss = 0.
      cs = 0.

      ccc = cos(freqs[i]*x)
      sss = sin(freqs[i]*x)
      xc = sum( y*ccc )
      xs = sum( y*sss )
      cc = sum( ccc**2 )
      ss = sum( sss**2 )
      cs = sum( ccc*sss )
      #for j in range(x.shape[0]):

          #c = cos(freqs[i] * x[j])
          #s = sin(freqs[i] * x[j])

          #xc += y[j] * c
          #xs += y[j] * s
          #cc += c * c
          #ss += s * s
          #cs += c * s

      tau = arctan(2 * cs / (cc - ss)) / (2 * freqs[i])
      c_tau = cos(freqs[i] * tau)
      s_tau = sin(freqs[i] * tau)
      c_tau2 = c_tau * c_tau
      s_tau2 = s_tau * s_tau
      cs_tau = 2 * c_tau * s_tau

      pgram[i] = 0.5 * (((c_tau * xc + s_tau * xs)**2 / \
          (c_tau2 * cc + cs_tau * cs + s_tau2 * ss)) + \
          ((c_tau * xs - s_tau * xc)**2 / \
          (c_tau2 * ss - cs_tau * cs + s_tau2 * cc)))

  return pgram

def computeSpectrum(x, y, freq) :
  # compute the spectrum
  pgram = lombscargle(x, y, 2*pi*freq)
  return sqrt(4*pgram/len(x))

def test() :
  # First define some input parameters for the signal:

  A = 2. # amplitude
  w = 1./10.0 # sampling frequency (Hz)
  phi = 0.5 * pi # phase
  nin = 1500
  nout = 1000
  frac_points = 0.9 # Fraction of points to discard

  # Randomly select a fraction of an array with timesteps:

  r = random.rand(nin)
  x = linspace(0.01, 100, nin)
  x_lossy = x[r >= frac_points]

  # Plot a sine wave for the selected times:

  y = A * sin(2*pi*w*x+phi)
  y_lossy = A * sin(2*pi*w*x_lossy+phi)

  # Define the array of frequencies for which to compute the periodogram:

  f = linspace(0.01, 1, nout)

  # Calculate Lomb-Scargle periodogram:

  #pg = lombscargle(x, y, 2*pi*f)
  pg = periodogram(x, y, f)
  pg_lossy = periodogram(x_lossy, y_lossy, f)

  # Now make a plot of the input data:
  import matplotlib.pyplot as plt
  import matplotlib.mlab
  pg_fft,f_fft = matplotlib.mlab.psd(y,NFFT=16*256,Fs=3./2*1./w) # only for evenly sampled data
  plt.subplot(2, 1, 1)
  plt.plot(x, y, 'r.-')
  plt.hold(True)
  plt.plot(x_lossy, y_lossy, 'bo')
  plt.xlabel('time (s)')

  # Then plot the normalized periodogram:

  ax2 = plt.subplot(2, 1, 2)
  #plt.plot(f, sqrt(4*(pg/len(x))))
  #plt.plot(f, pg/len(x))
  ax2.plot(1./f, pg,'r')
  plt.hold(True)
  ax2.plot(1./f, pg_lossy,'b')
  #ax2.plot(1./f_fft, sqrt(pg_fft),'g')
  ax2.set_xlim( [0,20] )
  ax2.legend(['100% data',str(int(round((1-frac_points)*100)))+'% data'])
  ax2.set_xlabel('period (s)')
  ax2.set_ylabel('amplitude')
  plt.show()
