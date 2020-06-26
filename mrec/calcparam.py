# -*- coding: utf-8 -*-

import xspec

def calcApecParam(
    rmf, arf, nH_1022pcm2, redshift, kT_keV, abun,
    emin_keV=0.5, emax_keV=2.0,
    lxemin_keV=0.1, lxemax_keV = 2.4,
    H0=72., L0=0.7
    ):
    """ Calculate parameters of an apec model using xspec
    
    Parameters:

    :rmf (str):
        Path to response matrix file
    
    :arf (str):
        Path to ancillary response file
    
    :nH_1022pcm2 (float):
        Total weighted galactic column density of Hydrogen in [10^22 atoms/cm**2]
    
    :redshift (float):
        Redshift of the cluster
    
    :kT_keV (float):
        Gas temperature in keV  
    
    :abun (float):
        Metalicity in solar units

    :emin_keV:
        Lower energy bound in which parameters are calculated
    
    :emax_keV:
        Upper  energy bound in which parameters are calculated
    
    :lxemin_keV:
        Lower energy bound in which luminosity is calculated
    
    :lxemax_keV:
        Upper  energy bound in which luminosity is calculated
    
    :H0 (float):
        Hubble constant in [km/s/Mpc]
    
    :L0 (float):
        Vacuum energy density


    Returns:
    
    Apec model count rate, flux and luminosity

    """

    # set xspec chattiness
    xspec.Xset.chatter = 0
  
    # define model
    m1 = xspec.Model("phabs*apec")

    # fake spectra
    fs = xspec.FakeitSettings(response=rmf, arf=arf, exposure=1e7)
    xspec.AllData.fakeit(nSpectra=1, settings=fs, noWrite=True)
  
    # spectrum  
    s1 = xspec.AllData(1)

    # set model parameters
    m1(1).values = nH_1022pcm2
    m1(2).values = kT_keV
    m1(3).values = abun
    m1(4).values = redshift
    m1(5).values = 1.


    # set abundance table to
    # Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197)
    xspec.Xset.abund = "angr"

    # set cosmology: cosmo <$H_0$><$q_0$><$\Lambda_0$>
    xspec.Xset.cosmo = "{0} .0 {1}".format(H0,L0)

    # set basic settings
    xspec.Plot.xAxis = "keV"
    s1.ignore("**-{0},{1}-**".format(emin_keV, emax_keV))

    # calculate model count rate
    modelCtr = s1.rate[3] # store predicted model rate in counts/s

    # calcualte flux in rest-frame band
    xspec.AllModels.calcFlux("{0} {1}".format(emin_keV, emax_keV))
    modelFlux = s1.flux[0] # store flux in erg/s/cm**2

    # calculate luminosity
    # since it's intrinsic luminosity, set nH to zero
    m1(1).values = 0.
    xspec.AllModels.calcLumin("{0} {1} {2}".format(lxemin_keV, lxemax_keV, redshift))
    modelLum = s1.lumin[0] * 1e44

    return modelCtr, modelFlux, modelLum