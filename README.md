# SIPK2
Systematics-Insensitive Periodograms for K2.

The paper associated with this code, "Systematics-insensitive periodic signal
search with K2" can be found at http://arxiv.org/abs/1505.07105.

To use this code, clone this repo and do

`from SIP import SIP`.

You'll need the Eigen light curves too.
You can download those from http://bbq.dfm.io/ketu/elcs.
You can also download any of DFM's ready-made raw K2 light curves from
http://bbq.dfm.io/ketu/lightcurves.

Usage:

`s2n, amp2s, w = SIP(x, y, basis, frequencies)`.

Returns the squared signal-to-noise-ratio, the squared amplitude and the
weights.
Check out the demo at
http://nbviewer.ipython.org/github/RuthAngus/SIPK2/blob/master/SIPdemo.ipynb.

# Licence

Copyright 2015 Ruth Angus.
Licensed under the terms of the MIT License (see LICENSE).
