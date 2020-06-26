# -*- coding: utf-8 -*-

"""

This is an example on how to run the mass reconstruction code.

"""

import sys

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

sys.path.append('../')
import mrec


def main():
    # define cosmology
    cosmo = FlatLambdaCDM(H0=72. * u.km / u.s / u.Mpc, Om0=0.3)

    # define path to files
    imgFile = 'EvtImg_0300_2200.fits[0]'
    expFile = 'ExpMap_0300_2200.fits[0]'
    bkgFile = None # The example is a simulation without background

    # aperture in arcmin that defines the cluster outskirts
    aperture = (2.5, 6.)
    
    # cluster properties
    ra = 57.3051
    dec = 39.184
    redshift = 0.192662

    # define conversion factor
    # energy range of the example is 0.3-2.2 keV
    lum_ctr = 1.46603383082e+45

    # initialize class
    mc = mrec.CalcMass(
        imgFile, expFile, bkgFile, aperture,
        ra, dec, redshift, lum_ctr, cosmo
        )

    # run computation
    outDir = './'
    mc.runFit(outDir, nwalkers=32, nstepsBurnin=100, nsteps=200)

    # true M500 from the simulation
    m500True = 5.759271e+14
    print('M500 (true): {0:.2e} * 10^14 Msolar'.format(m500True/1e14))
    
    # output
    #Mean acceptance fraction: 0.815
    #M500: (5.30e+00 -4.78e-01 +4.99e-01) * 10^14 Msolar
    #M500 (true): 5.76e+00 * 10^14 Msolar


if __name__ == "__main__":
    main()