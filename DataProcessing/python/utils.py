import os
import errno
import numpy as np
import uproot as up

class CacheManager:
    def __init__(self, name):
        self.name_ = name
        self.cache = up.ArrayCache("1 GB")

    def dump(self):
        obj = open(self.name_, 'wb')
        pickle.dump(self.cache, obj)
        obj.close()

    def load(self):
        try:
            obj = open(self.name_, 'rb')
            self.cache = pickle.load(obj, encoding='bytes')
            obj.close()
        except IOError:
            pass
        except EOFError:
            print('Perhaps the cache was not properly saved in a previous session?')
            raise
        return self.cache

def calculate_rect_side_for_plot_bins(bins, scale='linear', nlayers=28):
    nbins = len(bins)-1
    if scale == 'linear':
        factor = abs(bins[-1]-bins[0]) / nbins
        return factor * np.ones(nbins*nlayers)
    elif scale == 'log':
        centers = (bins[1:]+bins[:-1])/2
        distances = [(centers[1]-centers[0])/2]
        distances.extend( (centers[1:]-centers[:-1])/2 )
        return np.array(distances*nlayers)
    else:
        raise ValueError('height_for_plot_bins: Option not supported.')

def create_dir(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

def flatten_dataframe(df):
    flat_df = df.to_numpy().flatten()
    if isinstance(flat_df[0], (list,np.ndarray)):
        flat_df = [item for items in flat_df for item in items]
    return np.array(flat_df)

def get_mean_and_sigma(x, y=None):
    """calculate mean and std for 1D and 2D distributions"""
    if y is None:
        mean = np.mean(x)
        sigma = np.sqrt( np.mean(x - mean)**2 )
    else:
        mean_squared = np.sum(y*x**2)
        mean = np.sum(y*x)
        sumy = np.sum(y)
        if sumy == 0:
            mean = 0
            sigma = np.inf
        else:
            mean_squared /= sumy
            mean /= sumy
            sigma = np.sqrt( mean_squared - mean**2 )
    return mean, sigma

def get_layer_col(df, starts_with, ilayer=None):
    """Obtains list of columns that match a specific column name ending."""
    if ilayer is None:
        regex = '^'+starts_with
    else:
        regex = '^'+starts_with+'.*_layer'+str(ilayer)+'$'
    layer_col = df.columns.str.contains(regex, regex=True)
    if ilayer is not None:
        assert( single_true(layer_col) )
    return layer_col

def get_sigma_band(x, y=None):
    mean, sigma = get_mean_and_sigma(x, y)
    if mean==0 and sigma==np.inf:
        left = mean
        right = mean
    else:
        sigma /= np.sqrt(len(x))
        left, right = mean - sigma/2, mean + sigma/2
    return left, right

def height_for_plot_bins(bins, scale='linear', nlayers=28):
    nbins = len(bins)-1
    if scale == 'linear':
        factor = abs(bins[-1]-bins[0]) / nbins
        return factor * np.ones(nbins*nlayers)
    elif scale == 'log':
        centers = (bins[1:]+bins[:-1])/2
        distances = [(centers[1]-centers[0])/2]
        distances.extend( (centers[1:]-centers[:-1])/2 )
        return np.array(distances*nlayers)
    else:
        raise ValueError('height_for_plot_bins: Option not supported.')

def input_sanity_checks(flags, argv):
    """Checks the input arguments for obvious mistakes"""
    for elem in argv:
        if '--' in elem and elem[2:] not in flags.__dict__.keys():
            raise IOError('ERROR: You passed an undefined input argument!')
    if flags.showertype == 'had' and flags.datatype == 'sim_noproton':
        raise ValueError('There is no proton-free sample for hadronic showers.')

def single_true(iterable):
    """Checks if one and only one element of iterable is True"""
    i = iter(iterable)
    return any(i) and not any(i)

def peak_abciss(counts, edges):
    """Get x value of the maximum count"""
    m = np.max(counts)
    idx = np.where( counts == m )[0][0] #first element only in case there are multiple equal maxima
    x = (edges[idx] + edges[idx+1])/2
    return x

def print_input_data(files):
    print("Input data:")
    for x in files:
        if os.path.isfile(x):
            print(x)
        else:
            print(x + ' (not found)')

#import warnings
#def _warning(message, category = UserWarning, filename = '', lineno = -1):
#    print(message + ' (' + filename + ':' + str(lineno) + ')')
#warnings.showwarning = _warning
