# mrec

## Introduction

Python package to reconstruct masses of galaxy clusters by making use of the self-similar behavior of cluster outskirts (e.g., Zhang et al., 2007; Ghirardini et al., 2018a; Käfer et al., 2019), that is the universality of the shape of the X-ray emission. 

Illustrative flow chart to reconstruct the galaxy cluster mass using an MCMC posterior sampling technique with a Poisson likelihood to compare the measured aperture photon counts to the β-model predictions:

![Flow chart](./flowchart.pdf)

## Requirements

Requires:
* Numpy
* Scipy
* Astropy
* emcee
* PyXspec (optional)
