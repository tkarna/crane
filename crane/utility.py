"""
Frequently used functions.

Tuomas Karna 2015-08-21
"""
import datetime
import os
from crane import plt


def parseTimeStr(timeStr):
    """
    Parses time string to datetime object.

    Parameters
    ----------
    timeStr : string
        string to parse. Following formats are supported:
        %Y-%m-%d-%H-%M-%S'
        %Y-%m-%d-%H-%M
        %Y-%m-%d-%H
        %Y-%m-%d

    Returns
    -------
    dt : datetime
        datetime object
    """
    if len(timeStr.split('-')) == 3:
        dt = datetime.datetime.strptime(timeStr, '%Y-%m-%d')
    elif len(timeStr.split('-')) == 4:
        dt = datetime.datetime.strptime(timeStr, '%Y-%m-%d-%H')
    elif len(timeStr.split('-')) == 5:
        dt = datetime.datetime.strptime(timeStr, '%Y-%m-%d-%H-%M')
    elif len(timeStr.split('-')) == 6:
        dt = datetime.datetime.strptime(timeStr, '%Y-%m-%d-%H-%M-%S')
    else:
        raise Exception('Could not parse date string:', timeStr)
    return dt


def createDirectory(path):
    """
    Creates given directory if it does not exist already.

    Raises an exception if a file with the same name exists.
    """
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise Exception('file with same name exists', path)
    else:
        os.makedirs(path)
    return path


def saveFigure(
        path,
        filename,
        extensions,
        verbose=False,
        dpi=200,
        bbox_tight=False):
    """
    Saves currently open figure in the given format.
    If extensions is a list, figure is saved in multiple formats

    Parameters
    ----------
    path : string
        Directory where image is stored
    filename : string
        Output file name without an extensions
    extensions : string or a list of strings
        File extension(s) that define the image file format. Must be supported by matplotlib.
    verbose : bool
        If True prints file name on stdout
    dpi : int
        Output file dot density (default 200)
    bbox_tight : bool
        If set will autocrop the figure to fit content (even if larger than original figure).
    """
    kw = {}
    if bbox_tight:
        kw['bbox_inches'] = 'tight'
    if not isinstance(extensions, list):
        extensions = [extensions]
    for ext in extensions:
        f = os.path.join(path, filename + '.' + ext)
        if verbose:
            print 'saving to', f
        plt.savefig(f, dpi=dpi, **kw)
