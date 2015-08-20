"""
Various colormaps and related routines

Tuomas Karna 2014-06-04
"""
from crane.plotting import plotBase
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib

def cmap_discretize(cmap, N):
  """Return a discrete colormap from the continuous colormap cmap.

      cmap: colormap instance, eg. cm.jet.
      N: number of colors.

  Example
      x = resize(arange(100), (5,100))
      djet = cmap_discretize(cm.jet, 5)
      imshow(x, cmap=djet)
  """

  if type(cmap) == str:
    cmap = plt.get_cmap(cmap)
  colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
  colors_rgba = cmap(colors_i)
  indices = np.linspace(0, 1., N+1)
  cdict = {}
  for ki,key in enumerate(('red','green','blue')):
    cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
  # Return colormap object
  return colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))

    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]

    from matplotlib.colors import ListedColormap
    return ListedColormap(colors, cmap.name + "_grayscale")

def get_BlueRed():
  """
  Blue to red colormap with non-linear color transition
  Useful for plotting fluxes/velocity with emphasis on near zero values
  """
  cdict_BuRd = {'red':  ((0.0, 0.0, 0.0),
                    (0.28,0.1, 0.1),
                    (0.48, 0.9, 0.9),
                    (0.5, 1.0, 1.0),
                    (0.52, 1.0, 1.0),
                    (0.72,.95, .95),
                    (1.0, 0.4, 1.0)),

          'green': ((0.0, 0.0, 0.0),
                    (0.28,0.1, 0.1),
                    (0.48, .95, .95),
                    (0.5, 1.0, 1.0),
                    (0.52, 0.9, 0.9),
                    (0.72,0.1, 0.1),
                    (1.0, 0.0, 0.0)),

          'blue':  ((0.0, 0.0, 0.4),
                    (0.28,.95, .95),
                    (0.48, 1.0, 1.0),
                    (0.5, 1.0, 1.0),
                    (0.52, 0.9, 0.9),
                    (0.72,0.1, 0.1),
                    (1.0, 0.0, 0.0))
          }
  return colors.LinearSegmentedColormap('BlueRed', cdict_BuRd)

def get_BlueRed2():
  """
  A variant of the BlueRed colormap with discontinuity in the middle
  """
  cdict_BuRd2 = {'red':  ((0.0, 0.0, 0.0),
                    (0.25,0.1, 0.1),
                    (0.5, 0.9, 1.0),
                    (0.75,.95, .95),
                    (1.0, 0.4, 1.0)),

          'green': ((0.0, 0.0, 0.0),
                    (0.25,0.1, 0.1),
                    (0.5, .95, 0.9),
                    (0.75,0.1, 0.1),
                    (1.0, 0.0, 0.0)),

          'blue':  ((0.0, 0.0, 0.4),
                    (0.25,.95, .95),
                    (0.5, 1.0, 0.9),
                    (0.75,0.1, 0.1),
                    (1.0, 0.0, 0.0))
          }
  return colors.LinearSegmentedColormap('BlueRed2', cdict_BuRd2)

def get_DarkSpectral(N=256) :
  """
  A modified version of the Spectral_r colormap with decreased luminosity and increased saturation.
  """
  mapper = cm.ScalarMappable(cmap=plt.get_cmap('Spectral_r'))
  a = np.linspace(0,1,256)
  # get array of rgb colors (alpha == 1)
  rgb_colors = mapper.to_rgba(a)[:,:3]
  # convert to hsv
  hsv_colors = np.squeeze(colors.rgb_to_hsv(rgb_colors[:,None,:]))
  # replace v (value=luminosity)
  hsv_colors2 = hsv_colors.copy()
  s = hsv_colors2[:,1]
  v = hsv_colors2[:,2]
  v[v > 0.97] = 0.97 # darken yellows
  v *= np.linspace(0.8,1.15,len(v)) # darkned low end, lighten high end
  v[v > 1.00] = 1.0
  #s *= 1.1 # add saturation
  #s[s > 1.00] = 1.0
  #hsv_colors2[:,2] = np.ones_like(hsv_colors2[:,2])*0.8
  #hsv_colors2[:,1] = np.ones_like(hsv_colors2[:,2])*0.9
  # convert to rgb
  rgb_colors2 = np.squeeze(colors.hsv_to_rgb(hsv_colors2[:,None,:]))
  # make colormap
  newcmap = colors.LinearSegmentedColormap.from_list('DarkSpectral',rgb_colors2,N=N)
  return newcmap

def get_light_cubehelix(N=256, reverse=False):
  """A lighter variant of cubehelix colormap that doesn't reach black/white"""
  import imp
  try:
    imp.find_module('cubehelix')
    cubehelix_is_installed = True
  except ImportError:
    cubehelix_is_installed = False
  if cubehelix_is_installed:
    import cubehelix
    cmap = cubehelix.cmap(rot=-1.5, minLight=0.1, maxLight=0.92, gamma=0.8,
                          nlev=N,reverse=reverse)
    cmap.name = 'cubehelix_light'
    if reverse:
        cmap.name += '_r'
  else:
    print 'cubehelix colormap is not installed, reverting to Spectral_r'
    cmap = plt.get_cmap('Spectral_r')
  return cmap

def plotColormap(cmap,ax=None) :
  """Plots the colormap on given axis, if no axis given uses current axis."""
  if ax == None :
    ax = plt.gca()
  gradient = np.linspace(0, 1, 256)
  gradient = np.vstack((gradient, gradient))
  ax.imshow(gradient,cmap=cmap,aspect=10)
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_title(cmap.name)

def plotColormaps(cmap_list, fig=None) :
  """Plots all the given colormaps in the same plot."""
  N = len(cmap_list)
  if fig == None:
    fig=plt.figure(figsize=(8,1.2*N))
  for i in range(N) :
    ax = fig.add_subplot(N,1,i+1)
    plotColormap( cmap_list[i], ax )

def custom_div_cmap(numcolors=11, name='custom_div_cmap',
                    mincol='blue', midcol='white', maxcol='red'):
    """ Create a custom diverging colormap with three colors
    
    Default is blue to white to red with 11 colors.  Colors can be specified
    in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
    """

    from matplotlib.colors import LinearSegmentedColormap 
    
    cmap = LinearSegmentedColormap.from_list(name=name, 
                                             colors =[mincol, midcol, maxcol],
                                             N=numcolors)
    return cmap

def linear_Lab_cmap(colorList, values=None, N=11, name='custom_Lab_cmap',
                    minlightness=5, maxlightness=95, interpspace='rgb',
                    exp=1.0):
    """
    Creates a color map from a list of colors by forcing a monotonic lightness
    gradient across the colormap.

    The colormap is created by piecewise linear interpolation. If values are
    not given, the colors are distributed evenly across the colormap (i.e.
    between values 0.0 and 1.0).

    Colors are cast to Lab space:
    a,b are taken from the colors in colorList,
    L is linearly interpolated between minlightness, maxlightness

    Parameters
    ----------

    colorList : list of colors
        list of colors in any format colorConverter.to_rgb() supports
    values : list of floats (optional)
        list of increasing values between 0.0 and 1.0 corresponding to
        the colors in the list
    N : int
        number of colors in the final map
    nane : str
        name of the colormap
    minlightness, maxlightness : float
        min/max value of lightness in Lab space, between 0 and 100
    interpspace : "rgb" | "lab"
        color space where the colors will be interpolated. "rgb" is normally
        better, "lab" may produce better transition for simple colormaps.
        Monotonic lightness can be guaranteed only with "lab" interpolation.
    exp : float
        Instead of varying lightness linearly, can apply power law:
        L = value**exp

    Returns
    -------

    cmap, cmap_r
        colormap and reversed colormap

    Tuomas Karna 2014-12-30
    """
    from skimage import color as skicolor
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap 
    from scipy.interpolate import interp1d

    cList = list(colorList)
    if values is None:
        values = np.linspace(0,1,len(cList))
    if values[0] > 0.0:
        values = np.hstack(([0], values))
        cList.insert(0, cList[0])
    if values[-1] < 1.0:
        values = np.hstack((values, [1.0]))
        cList.insert(len(cList), cList[-1])
    
    cList = [list(colors.colorConverter.to_rgb(c)) for c in cList]
    nPoints = len(cList)
    cList_Lab = skicolor.rgb2lab([cList])[0]
    cList_Lab = np.array(cList_Lab)

    if interpspace == 'rgb':
        # replace L channel by lin. gradient and create LinearSegmentedColormap
        L = interp1d(np.linspace(0,1,2), [minlightness, maxlightness])(values)
        # Steven's power law
        L = 100*((L/100)**exp)
        colors_Lab = cList_Lab.copy()
        colors_Lab[:, 0] = L
        colors_rgb = skicolor.lab2rgb([colors_Lab])[0]
        # crop bad RGB values
        colors_rgb[colors_rgb < 0.0] = 0.0
        colors_rgb[colors_rgb > 1.0] = 1.0
        cmap = LinearSegmentedColormap.from_list(name=name,
                                                 colors=colors_rgb,
                                                 N=N)
        cmap_r = LinearSegmentedColormap.from_list(name=name+'_r',
                                                 colors=colors_rgb[::-1],
                                                 N=N)
    else:
        # interpolate whole color array in Lab space
        x = np.linspace(0,1,N)
        L = interp1d(np.linspace(0,1,2), [minlightness, maxlightness])(x)
        # Steven's power law
        L = 100*((L/100)**exp)
        a = interp1d(values, cList_Lab[:,1])(x)
        b = interp1d(values, cList_Lab[:,2])(x)
        colors_Lab = np.vstack((L,a,b)).T
        colors_rgb = skicolor.lab2rgb([colors_Lab])[0]
        # crop bad RGB values
        colors_rgb[colors_rgb < 0.0] = 0.0
        colors_rgb[colors_rgb > 1.0] = 1.0
        cmap = ListedColormap(colors_rgb, name=name, N=None)
        cmap_r = ListedColormap(colors_rgb[::-1], name=name+'_r', N=None)

    return cmap, cmap_r

def linear_cmap(colorList, values=None, N=11, name='custom_Lab_cmap',
                interpspace='rgb'):
    """
    Creates a color map from a list of colors with linear interpolation.

    The colormap is created by piecewise linear interpolation. If values are
    not given, the colors are distributed evenly across the colormap (i.e.
    between values 0.0 and 1.0).

    Parameters
    ----------

    colorList : list of colors
        list of colors in any format colorConverter.to_rgb() supports
    values : list of floats (optional)
        list of increasing values between 0.0 and 1.0 corresponding to
        the colors in the list
    N : int
        number of colors in the final map
    nane : str
        name of the colormap
    interpspace : "rgb" | "lab"
        color space where the colors will be interpolated. "rgb" is normally
        better, "lab" may produce better transition for simple colormaps.
        Monotonic lightness can be guaranteed only with "lab" interpolation.

    Returns
    -------

    cmap, cmap_r
        colormap and reversed colormap

    Tuomas Karna 2014-12-30
    """
    from skimage import color as skicolor
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap 
    from scipy.interpolate import interp1d

    cList = list(colorList)
    if values is None:
        values = np.linspace(0,1,len(cList))
    if values[0] > 0.0:
        values = np.hstack(([0], values))
        cList.insert(0, cList[0])
    if values[-1] < 1.0:
        values = np.hstack((values, [1.0]))
        cList.insert(len(cList), cList[-1])

    cList = [list(colors.colorConverter.to_rgb(c)) for c in cList]
    nPoints = len(cList)
    cList_Lab = skicolor.rgb2lab([cList])[0]
    cList_Lab = np.array(cList_Lab)

    if interpspace == 'rgb':
        colors_Lab = cList_Lab.copy()
        colors_rgb = skicolor.lab2rgb([colors_Lab])[0]
        # crop bad RGB values
        colors_rgb[colors_rgb < 0.0] = 0.0
        colors_rgb[colors_rgb > 1.0] = 1.0
        cmap = LinearSegmentedColormap.from_list(name=name,
                                                 colors=colors_rgb,
                                                 N=N)
        cmap_r = LinearSegmentedColormap.from_list(name=name+'_r',
                                                 colors=colors_rgb[::-1],
                                                 N=N)
    else:
        # interpolate whole color array in Lab space
        x = np.linspace(0,1,N)
        L = interp1d(values, cList_Lab[:,0])(x)
        a = interp1d(values, cList_Lab[:,1])(x)
        b = interp1d(values, cList_Lab[:,2])(x)
        colors_Lab = np.vstack((L,a,b)).T
        colors_rgb = skicolor.lab2rgb([colors_Lab])[0]
        # crop bad RGB values
        colors_rgb[colors_rgb < 0.0] = 0.0
        colors_rgb[colors_rgb > 1.0] = 1.0
        cmap = ListedColormap(colors_rgb, name=name, N=None)
        cmap_r = ListedColormap(colors_rgb[::-1], name=name+'_r', N=None)

    return cmap, cmap_r

cm.register_cmap(cmap=get_light_cubehelix())
cm.register_cmap(cmap=get_light_cubehelix(reverse=True))


def reverse_colormap(cmap):
    assert isinstance(cmap, matplotlib.colors.LinearSegmentedColormap), \
        'colormap must be LinearSegmentedColormap'
    rgb_r = []
    for channel in ['red', 'green', 'blue']:
        data = []
        for t in cmap._segmentdata[channel]:
            data.append((1.0 - t[0], t[1], t[2]))
        rgb_r.append(sorted(data))
    segmentdata = dict(zip(['red', 'green', 'blue'], rgb_r))
    name = cmap.name
    cmap_r = matplotlib.colors.LinearSegmentedColormap(name + '_r', segmentdata)
    return cmap_r