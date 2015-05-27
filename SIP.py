import numpy as np

# loop over eval_freq to compute a periodogram
def SIP(x, y, basis, fs):
    """
    # `K2pgram`

    Calls eval_freq on x and y for each frequency in fs.
    fs is frequency, *not* angular freq.
    Returns the sum of squared amplitudes and the signal to noise.
    """
    # construct arrays
    AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
    ATA = np.dot(AT, AT.T)

    amps2 = np.zeros_like(fs)
    s2n = np.zeros_like(fs)
    for i, f in enumerate(fs):
        amps2[i], s2n[i], w = eval_freq(x, y, f, AT, ATA)
    return s2n, amps2, w

# calculate periodogram by just updating sin and cos parts of the  matrices
def eval_freq(x, y, f, AT, ATA, compute_trends=False):
    """
    # `eval_freq`

    Computes the sum of squared amplitudes for a sine and cosine function
    with frequency, f.
    Linear least squares fit to the data.
    """
    arg = 2*np.pi*f*x
    AT[-2, :] = np.sin(arg)
    AT[-1, :] = np.cos(arg)

    # AT.shape = (153, nt)
    # shape: (151, nt) * (nt, 2) -> (151, 2)
    v = np.dot(AT[:-2, :], AT[-2:, :].T)
    ATA[:-2, -2:] = v
    ATA[-2:, :-2] = v.T

    # AT[-2:, :].shape = (2, nt)
    # (2, nt), (nt, 2)
    ATA[-2:, -2:] = np.dot(AT[-2:, :], AT[-2:, :].T)
    w = np.linalg.solve(ATA, np.dot(AT, y))

    S = np.linalg.inv(ATA)[-2:, -2:]
    s2n = np.dot(w[-2:], np.linalg.solve(S, w[-2:]))

    if compute_trends:
        trends = np.dot(w[:-2], AT[:-2])
        return np.sum(w[-2:]**2), s2n, trends
    return s2n, np.sum(w[-2:]**2), w

# compute trends at a certain frequency.
def eval_1_freq(x, y, f, compute_trends=True):
    """
    # `eval_1_freq`

    Computes the sum of squared amplitudes for a sine and cosine function
    with frequency, f.
    Doesn't use precomputed matrices so is slow compared to eval_freq.
    This version is just for computing conditional light curves.
    It will provide the same output as eval_freq if you give it the
    pre-constructed matrices.
    """

    # construct arrays
    AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
    ATA = np.dot(AT, AT.T)

    arg = 2*np.pi*f*x
    AT[-2, :] = np.sin(arg)
    AT[-1, :] = np.cos(arg)

    # AT.shape = (153, nt)
    # shape: (151, nt) * (nt, 2) -> (151, 2)
    v = np.dot(AT[:-2, :], AT[-2:, :].T)
    ATA[:-2, -2:] = v
    ATA[-2:, :-2] = v.T

    # AT[-2:, :].shape = (2, nt)
    # (2, nt), (nt, 2)
    ATA[-2:, -2:] = np.dot(AT[-2:, :], AT[-2:, :].T)
    w = np.linalg.solve(ATA, np.dot(AT, y))

    S = np.linalg.inv(ATA)[-2:, -2:]
    s2n = np.dot(w[-2:], np.linalg.solve(S, w[-2:]))

    if compute_trends:
        trends = np.dot(w[:-2], AT[:-2])
        return np.sum(w[-2:]**2), s2n, trends
    return s2n, np.sum(w[-2:]**2), w
