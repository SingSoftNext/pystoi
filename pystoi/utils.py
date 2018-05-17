import numpy as np
import scipy

EPS = np.finfo("float").eps

def thirdoct(fs, nfft, num_bands, min_freq):
    """ Returns the 1/3 octave band matrix and its center frequencies
    # Arguments :
        fs : sampling rate
        nfft : FFT size
        num_bands : number of 1/3 octave bands
        min_freq : center frequency of the lowest 1/3 octave band
    # Returns :
        obm : Octave Band Matrix
        cf : center frequencies
    """
    f = np.linspace(0, fs, nfft + 1)
    f = f[:int(nfft/2) + 1]
    k = np.array(range(num_bands)).astype(float)
    cf = np.power(2. ** (1. / 3), k) * min_freq
    freq_low = min_freq * np.power(2., (2 * k -1 ) / 6)
    freq_high = min_freq * np.power(2., (2 * k + 1) / 6)
    obm = np.zeros((num_bands, len(f))) # a verifier

    for i in range(len(cf)):
        # Match 1/3 oct band freq with fft frequency bin
        f_bin = np.argmin(np.square(f - freq_low[i]))
        freq_low[i] = f[f_bin]
        fl_ii = f_bin
        f_bin = np.argmin(np.square(f - freq_high[i]))
        freq_high[i] = f[f_bin]
        fh_ii = f_bin
        # Assign to the octave band matrix
        obm[i, fl_ii:fh_ii] = 1
    return obm, cf


def stft(x, win_size, fft_size, overlap=4):
    """ Short-time Fourier transform for real 1-D inputs
    # Arguments
        x : 1D array, the waveform
        win_size : integer, the size of the window and the signal frames
        fft_size : integer, the size of the fft in samples (zero-padding or not)
        overlap: integer, number of steps to make in fftsize
    # Returns
        stft_out : 2D complex array, the STFT of x.
    """
    hop = int(win_size / overlap)
    w = scipy.hanning(win_size + 2)[1: -1]  # = matlab.hanning(win_size)
    stft_out = np.array([np.fft.rfft(w * x[i:i + win_size], n=fft_size)
                        for i in range(0, len(x) - win_size, hop)])
    return stft_out


def remove_silent_frames(x, y, dyn_range, framelen, hop):
    """ Remove silent frames of x and y based on x
    A frame is excluded if its energy is lower than max(energy) - dyn_range
    The frame exclusion is based solely on x, the clean speech signal

    # Arguments :
        x : array, original speech wav file
        y : array, denoised speech wav file
        dyn_range : Energy range to determine which frame is silent
        framelen : Window size for energy evaluation
        hop : Hop size for energy evaluation

    # Returns :
        x without the silent frames
        y without the silent frames (aligned to x)
    """
    # Compute Mask
    #!!!
    x = np.zeros(16)
    framelen = 8
    hop = 4
    noverlap = framelen - hop

    # Pad x to integer number of bins and for endpoint reconstruction
    pad = ...
    #!!!
    w = scipy.hanning(framelen + 2)[1:-1]
    #!!! Testing reconstruction
    from scipy import signal
    w = np.sqrt(signal.hann(framelen, False))
    w = signal.hann(framelen, False)
    print(f'COLA: {signal.check_COLA(w, framelen, noverlap)}')
    x = np.random.uniform(-1, 1, len(x) + 2 * framelen)
    print(f'x old: {len(x)}')
    nbins = 1 + (len(x) - framelen) // hop
    x = x[:(nbins - 1) * hop + framelen]  # eliminate partial bins
    print(f'x new: {len(x)}')
    x[:framelen] = 0
    x[-framelen:] = 0
    print(f'#bins: {1 + (len(x) - framelen) // hop}')
    y = x

    # Calculate normalization so overlap-add works
    ola_denom = np.zeros_like(x)
    for i in range(nbins):
        ola_denom[i*hop:i*hop+framelen] += w
    #!!!
    # mask = np.array([20 * np.log10(np.linalg.norm(w * x[i:i + framelen]) + EPS)
    #                  for i in range(0, len(x) - framelen, hop)])
    mask = np.array([20 * np.log10(np.linalg.norm(w * x[i:i + framelen]) + EPS)
                     for i in range(0, len(x) - framelen + 1, hop)])
    mask += dyn_range - np.max(mask)
    mask = mask > 0
    #!!!
    mask[:] = True
    print(f'len(mask): {len(mask)}')
    print(f'#True:     {np.sum(mask)}')
    #!!!
    # Remove silent frames
    count = 0
    x_sil = np.zeros(x.shape)
    y_sil = np.zeros(y.shape)
    for i in range(mask.shape[0]):
        if mask[i]:
            sil_ind = range(count * hop, count * hop + framelen)
            ind = range(i * hop, i * hop + framelen)
            x_sil[sil_ind] += x[ind] * w
            y_sil[sil_ind] += y[ind] * w
            count += 1
    # Cut unused length
    x_sil = x_sil[:sil_ind[-1] + 1]
    y_sil = y_sil[:sil_ind[-1] + 1]
    #!!!
    x = x[framelen:-framelen]
    x_sil = x_sil[framelen:-framelen]
    ola_denom = ola_denom[framelen:-framelen]
    print(f'ola_demom avg: {np.mean(ola_denom):0.2f}')
    print(f'ola_denom std: {np.std(ola_denom):0.2f}')
    x_sil /= ola_denom
    print(f'Close?: {np.allclose(x, x_sil)}')
    print(f'delta max: {np.max(np.abs(x - x_sil)):0.2f}')
    print(f'arg max: {np.argmax(np.abs(x - x_sil))}')
    exit()
    #!!!
    return x_sil, y_sil


def corr(x, y):
    """ Returns correlation coefficient between x and y (1-D)"""
    new_x = x - np.mean(x)
    new_x /= np.sqrt(np.sum(np.square(new_x)))
    new_y = y - np.mean(y)
    new_y /= np.sqrt(np.sum(np.square(new_y)))
    rho = np.sum(new_x * new_y)
    return rho


def vect_two_norm(x, axis=-1):
    """ Returns a vectors of norms of the rows of a matrix """
    return np.sum(np.square(x), axis=axis)


def row_col_normalize(x):
    """ Row and column mean and variance normalize a 2D matrix """
    # Row mean and variance normalization
    x_normed = x + EPS * np.random.standard_normal(x.shape)
    x_normed -= np.mean(x_normed, axis=-1, keepdims=True)
    x_normed = np.matmul(np.diag(1. / np.sqrt(vect_two_norm(x_normed))), x_normed)
    # Column mean and variance normalization
    x_normed += + EPS * np.random.standard_normal(x_normed.shape)
    x_normed -= np.mean(x_normed, axis=0, keepdims=True)
    x_normed = np.matmul(x_normed, np.diag(1. / np.sqrt(vect_two_norm(x_normed, axis=0))))
    return x_normed
