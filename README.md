# SMAP
This repository contains Julia and Python notebooks (jupyter) for some of the examples in my book "Stochastic Methods in Asset Pricing" (SMAP), MIT Press 2017. Among other things, the method for computing the price of American call options and the construction of the early exercise premium curve in the Black-Scholes-Merton framework, as described in section 18.4 in SMAP, is included in both Python and Julia versions. This approach is very different from the classical Monte Carlo, binomial lattice, or PDE methods. I hope that someone could find a way to speed up the code.

Some more mundane task, such as storing an uploading market data, simulating multivariate normal samples with a given covariance matrix, building historams, and such, are implemented as well. For some reason I was not able to get the normalization in the 'histogram' function in Julia to do what is expected, so I wrote my one histogram function, which I call 'hstgram'. This matter requires a follow-up.   


