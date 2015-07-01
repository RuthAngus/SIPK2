import numpy as np
import matplotlib.pyplot as plt
import h5py
import fitsio
from SIP import SIP, eval_freq

def EPICSIP(EPIC, C, periods, nELC=150, plot=False):
    '''
    Given an EPIC ID, campaign number and period array,
    produce a sip.
    '''
    fname = "ktwo%s-c%s_lpd-lc.fits" % (str(int(EPIC)), str(int(C)).zfill(2))
    print fname, "found"

    # load the data
    data = fitsio.read(fname)
    aps = fitsio.read(fname, 2)
    y = data["flux"][:, np.argmin(aps["cdpp6"])]
    x = data["time"]
    q = data["quality"]
    m = np.isfinite(y) * np.isfinite(x) * (q==0)
    y, x = y[m], x[m]
    y = y / np.median(y) - 1

    nELC = 150
    with h5py.File("c%s.h5" % C, "r") as f:
        basis = f["basis"][:nELC, m]

    if plot:
        plt.clf()
        plt.plot(x, y, "k")
        plt.xlabel("BJD - 2454833 (days)")
        plt.ylabel("Normalised Flux")
        plt.savefig("%s_lc" % EPIC)

    freqs = 1./periods
    s2n, amp2s, w = SIP(x, y, basis, freqs)

    if plot:
        plt.clf()
        plt.plot(periods, amp2s/max(amp2s), "r", label="amp2s")
        plt.plot(periods, s2n/max(s2n), "k", label="s2n")
        plt.legend()
        plt.xlabel("Period (days)")
        plt.ylabel("Relative (S/N)^2")
        plt.savefig("%s_sip" % EPIC)

    return x, y, s2n, amp2s, w
