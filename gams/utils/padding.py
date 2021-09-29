import numpy as np


def derivative_mirror(grid_vals, direction, pads, mirror='inverse'):
    if mirror.lower() == 'inverse':
        sign = 1
    else:
        sign = -1
    nx, ny = grid_vals.shape
    pad_left, pad_right, pad_bot, pad_top = pads
    # if ny > nx:
    if direction == 'lr':
        deriv = np.diff(grid_vals, axis=0)
        pad_vals = np.pad(grid_vals, [[pad_left, pad_right], [0, 0]], mode='constant')
        # deriv = np.diff(pad_vals, axis=0)
        # for ii in range(pad_left):
        #     pad_vals[pad_left-1-ii,:] = pad_vals[pad_left-ii, :] + deriv[pad_left+ii+1, :]
        # for ii in range(pad_right-1):
        #     pad_vals[pad_left+nx+ii+1,:] = pad_vals[pad_left+nx+ii, :] + deriv[pad_left+nx-ii-2, :]
        for ii in range(pad_left):
            pad_vals[pad_left-1-ii,:] = pad_vals[pad_left-ii, :] - sign*deriv[ii, :]
        for ii in range(pad_right):
            pad_vals[pad_left+nx+ii,:] = pad_vals[pad_left+nx+ii-1, :] + sign*deriv[nx-ii-2, :]
    elif direction == 'ud':
        deriv = np.diff(grid_vals, axis=1)
        pad_vals = np.pad(grid_vals, [[0, 0], [pad_bot, pad_top]])
        # deriv = np.diff(pad_vals, axis=1)
        for ii in range(pad_bot):
            # pad_vals[:,pad_bot-1-ii] = pad_vals[:, pad_bot-ii] + deriv[:, pad_bot+ii+1]
            pad_vals[:,pad_bot-1-ii] = pad_vals[:, pad_bot-ii] - sign*deriv[:, ii]
        for ii in range(pad_top):
            # pad_vals[:, pad_bot+ny+ii] = pad_vals[:, pad_bot+ny+ii-1] + deriv[:,pad_bot+ny-ii-1]
            pad_vals[:, pad_bot+ny+ii] = pad_vals[:, pad_bot+ny+ii-1] + sign*deriv[:,ny-ii-2]
    return pad_vals


def u2v(u: np.ndarray) -> np.ndarray:
    '''Converts the image `u` into the image `v`
    Parameters
    ----------
    u : np.ndarray
        [M, N] image
    Returns
    -------
    v : np.ndarray
        [M, N] image, zeroed expect for the outermost rows and cols
    '''
    v = np.zeros(u.shape, dtype=np.float64)

    v[0, :] = np.subtract(u[-1, :], u[0,  :], dtype=np.float64)
    v[-1,:] = np.subtract(u[0,  :], u[-1, :], dtype=np.float64)

    v[:,  0] += np.subtract(u[:, -1], u[:,  0], dtype=np.float64)
    v[:, -1] += np.subtract(u[:,  0], u[:, -1], dtype=np.float64)
    return v


def v2s(v_hat: np.ndarray) -> np.ndarray:
    '''Computes the maximally smooth component of `u`, `s` from `v`
    s[q, r] = v[q, r] / (2*np.cos( (2*np.pi*q)/M )
        + 2*np.cos( (2*np.pi*r)/N ) - 4)
    Parameters
    ----------
    v_hat : np.ndarray
        [M, N] DFT of v
    '''
    M, N = v_hat.shape

    q = np.arange(M).reshape(M, 1).astype(v_hat.dtype)
    r = np.arange(N).reshape(1, N).astype(v_hat.dtype)

    den = (2*np.cos( np.divide((2*np.pi*q), M) ) \
         + 2*np.cos( np.divide((2*np.pi*r), N) ) - 4)
    # s = np.divide(v_hat, den, out=np.zeros_like(v_hat), where=den!=0.0)
    s = np.zeros(v_hat.shape, dtype=np.complex128)
    idx = den != 0
    s[idx] = v_hat[idx] / den[idx]
    s[0, 0] = 0
    # s = den
    return s


def periodic_smooth_decomp(I: np.ndarray) -> (np.ndarray, np.ndarray):
    '''Performs periodic-smooth image decomposition
    Parameters
    ----------
    I : np.ndarray
        [M, N] image. will be coerced to a float.
    Returns
    -------
    P : np.ndarray
        [M, N] image, float. periodic portion.
    S : np.ndarray
        [M, N] image, float. smooth portion.
    '''
    u = I.astype(np.float64)
    v = u2v(u)
    v_fft = np.fft.fftn(v)
    s = v2s(v_fft)
    s_i = np.fft.ifftn(s)
    s_f = np.real(s_i)
    p = u - s_f # u = p + s
    return p, s_f
    # u = I
    # return u, u


def zeros(grid_vals, padding):
    return np.pad(grid_vals, padding, mode='constant', constant_values=0)


def wrap(grid_vals, padding):
    return np.pad(grid_vals, padding, mode='wrap')


def reflect(grid_vals, padding, inverse=False):
    pad_vals = np.pad(grid_vals, padding, mode='reflect')
    if inverse:
        pad_left, pad_right = padding[0]
        pad_bot, pad_top = padding[1]
        nx, ny = grid_vals.shape
        idx_right = pad_left + nx
        idx_top = pad_bot + ny
        # Left fold
        if pad_left:
            val = 2*pad_vals[pad_left+2,pad_bot:pad_bot+ny]
            pad_vals[:pad_left,pad_bot:pad_bot+ny] = -1*pad_vals[:pad_left,pad_bot:idx_top] + val
        # right fold
        if pad_right:
            val = 2*pad_vals[idx_right,pad_bot:pad_bot+ny]
            pad_vals[idx_right:,pad_bot:idx_top] = -1*pad_vals[idx_right:,pad_bot:idx_top] + val
        # top fold
        if pad_top:
            val = 2*pad_vals[pad_left:idx_right, idx_top]
            pad_vals[pad_left:idx_right,pad_bot+ny:] = -1*pad_vals[pad_left:idx_right,idx_top:] + val[:, np.newaxis]
        # bottom fold
        if pad_bot:
            val = 2*pad_vals[pad_left:pad_left+nx, pad_bot+2]
            pad_vals[pad_left:pad_left+nx, :pad_bot] = -1*pad_vals[pad_left:pad_left+nx, :pad_bot] + val[:, np.newaxis]
    return pad_vals


def psd(grid_vals, padding=None):
    if padding:
        # padding = [[0, 0], [0, 0]]
        direction = 'ud' if grid_vals.shape[0] > grid_vals.shape[1] else 'lr'
        pad_left, pad_right = padding[0]
        pad_bot, pad_top = padding[1]
        grid_vals = derivative_mirror(grid_vals, direction, [pad_left, pad_right, pad_bot, pad_top])
    # pad_vals = periodic_smooth_decomp(pad_vals)[0]
    pad_vals = periodic_smooth_decomp(grid_vals)[0]
    # pad_vals = grid_vals
    return pad_vals