# Epidemic-Transport-Model

In this repository, we publish the code that has been used to perform the numerical simulations in

C. Kuehn & J. Mölter. "The influence of a transport process on the epidemic threshold". *J. Math. Biol.* **85**:62 (2022). DOI: [10.1007/s00285-022-01810-7](https://doi.org/10.1007/s00285-022-01810-7). arXiv:[2112.04951 [nlin.AO]](https://doi.org/10.48550/arXiv.2112.04951).

## Overview

We consider a model of an epidemic in the presence of a transport process. In this model, the epidemic dynamics are supported by a two-layer multiplex network structure. The bottom layer is a static network representing the social relations underlying a population whereas the top layer is dynamically changing its topology subject to a transport process taking place on a separate transport network. More precisely, the individuals of the population perform a random walk on the transport network and whether two of them occupy the same site they become linked in the dynamic layer of the epidemic network until one of them leaves the site.

## Code

We provide both C++ and Python implementations to simulate this model. As far as the event-based stochastic simulation of the epidemic spreading goes, these are essentially equivalent. However, the Python implementation in addition has the capability to compute the corresponding mean-field solutions.

Note: Besides the standard library, the C++ implementation builds on the [igraph](https://igraph.org/) as well as the [Armadillo](http://arma.sourceforge.net/) libraries.

## References

C. Kuehn & J. Mölter. "The influence of a transport process on the epidemic threshold". *J. Math. Biol.* **85**:62 (2022). DOI: [10.1007/s00285-022-01810-7](https://doi.org/10.1007/s00285-022-01810-7). arXiv:[2112.04951 [nlin.AO]](https://doi.org/10.48550/arXiv.2112.04951).
