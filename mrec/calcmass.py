# -*- coding: utf-8 -*-

import numpy as np
import os
import emcee
import math

from scipy.integrate import quad

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord



class CalcMass:
    """
    Class to reconstruct galaxy cluster masses
    based on the surface brightness in the outskirts
    """

    def __init__(
        self, imgFile, expFile, bkgFile, aperture,
        ra, dec, redshift, lum_ctr, cosmo
        ):

        """
        Parameters:

        :imgFile (str):
            Path to file that stores the photon image
        
        :expFile (str):
            Path to file that stores the exposure map (optional)
        
        :bkgFile (str):
            Path to file that stores the background image (optional)

        :aperture (None or tuple of floats):
            Inner and outer aperture bounds in arcmin, e.g. (2., 5.)
        
        :ra (float):
            Right ascension of the cluster in degrees

        :dec (float):
            Declination of the cluster in degrees

        :redshift (float):
            Redshift of the cluster

        :lum_ctr (float):
            Conversion factor to transform luminosity into count rates
        
        :cosmo:
            Instance of the astropy.cosmology sub-package
        """ 

        self.imgFile =imgFile
        self.expFile = expFile
        self.bkgFile = bkgFile
        self.aperture = aperture

        self.ra = ra
        self.dec = dec
        self.redshift = redshift
        self.lum_ctr = lum_ctr
        
        self.cosmo = cosmo
        self.ez = cosmo.efunc(redshift)
        self.rhoCritz = cosmo.critical_density(redshift).value # g / cm3
        self.sf_kpcparcmin = cosmo.kpc_proper_per_arcmin(redshift).value

        # define physical constants and conversion factors
        self.msun = 1.98892e33 # g/Msolar
        self.kpc_cm = 3.085678e21

        self.prepFiles()


    def loadFitsImage(self, finput, attrStr):
        """ Function to check the validity of the input fits file
        
        Input:

        :finput (str): Path to image file, e.g. "./examples/evt.fits[0]"
        :attrStr (str): Unique name for attributes, e.g. "img"
        
        Returns:

        header and image instances:
        accessible via self.{attStr}Hdr, self.{attrStr}Data
        """

        if finput[-1] == ']':
            ext = int(finput.split('[')[-1].strip(']'))
            finput = finput.split('[')[0]
        else:
            #print('No extension number provided. It is set to 0')
            ext = 0

        if not os.path.isfile(finput):
            raise FileNotFoundError('Input fits file not found: {0}'.format(finput))

        with fits.open(finput) as hdul:
            if len(hdul)-1 < ext:
                raise FileNotFoundError('Provided extension number {0} not found'.format(ext))
    
            hdr  = hdul[ext].header

            if ext == 0 and hdr['NAXIS'] == 2:
                data = hdul[ext].data
            elif ext != 0 and hdr['XTENSION'] == 'IMAGE':
                data = hdul[ext].data
            else:
                raise FileNotFoundError('No IMAGE found in extension {0}'.format(ext))
        
        setattr(self, '{0}Hdr'.format(attrStr), hdr)
        setattr(self, '{0}Data'.format(attrStr), data)


    def prepFiles(self):  
        # load cluster image, exposure, and background map
        self.loadFitsImage(self.imgFile, 'img')
        
        if self.expFile is not None:
            self.loadFitsImage(self.expFile, 'exp')
        else:
            self.expHdr, self.expData = None, np.ones(np.shape(self.imgData))
            
        if self.bkgFile is not None:
            self.loadFitsImage(self.bkgFile, 'bkg')
        else:
            self.bkgHdr, self.bkgData = None, np.zeros(np.shape(self.imgData))

        # image WCS and cluster position
        self.imgWcs = WCS(self.imgFile.split('[')[0])
        self.cluSky = SkyCoord(self.ra, self.dec, unit='deg', frame='fk5')
        self.calcPixelSeparation()

        # calcualte mask if flag is not set to self calibration
        if self.aperture != None:
            self.regMask = np.logical_and(
                self.pixSep >= self.aperture[0],
                self.pixSep <= self.aperture[1]
                )

        # get pixel size
        try:
            self.pixSize = self.imgHdr['CDELT2']
            self.pixUnit = self.imgHdr['CUNIT2']       
        except KeyError:
            self.pixSize = 4./60.
            self.pixUnit = 'arcmin'
            print('Pixel size/unit not in image header.')
            print('Assume: {0} {1}'.format(self.pixSize, self.pixUnit))


    def calcPixelSeparation(self):
        """ Calcualte a map storing pixel separations in arcmin """
        imgShape = np.shape(self.imgData)
    
        iRaDec = np.array([
            [iRa,iDec] for iRa in range(imgShape[1]) for iDec in range(imgShape[0])
            ])
        
        raDec = self.imgWcs.all_pix2world(iRaDec[:,0], iRaDec[:,1], 1, ra_dec_order=True)
        pixSky = SkyCoord(raDec[0], raDec[1], unit='deg', frame='fk5')
        pixSep = self.cluSky.separation(pixSky).arcminute
        self.pixSep = pixSep.reshape(np.shape(self.imgData))


    # single beta model w/o background
    # normalization is set to 1.
    def singleBeta(self, x, beta, rc):
        return ( 1. + (x/rc)**2. ) ** ( -3. * beta + 0.5)

    def crIntegrandSingleBeta(self, x, beta, rc):
        return self.singleBeta(x, beta, rc) * x

    def intSingleBetaModel(self, apIn, apOut, beta, rc):
        return 2.*np.pi*quad(self.crIntegrandSingleBeta, apIn, apOut, args=(beta, rc))[0]


    def calcLx(self, m5):
        """ Use kettula et al 2015 M500-Lx relation
        Lx measured in 0.1-2.4 keV band, inside 0.1-1 r500
        """
        alpha_ml = 0.7
        log10N_ml = 0.15

        lx = 1e44*self.ez * 10**( (np.log10(m5 *self.ez/3e14) - log10N_ml) / alpha_ml)
        return lx # erg/s


    def calcR5(self, m5):
        r5 = np.power(np.array(m5)*self.msun/4*3/np.pi/500/self.rhoCritz, 1./3) # cm
        return r5 / self.kpc_cm / self.sf_kpcparcmin # arcmin
  

    def calcPercentiles(self, a,b,c):
        """a,b,c: 15.87, 50., 84.13 percentiles
        """
        a,b,c = np.array(a),np.array(b),np.array(c)
        return np.array([b, b-a, c-b])


    def lnlike(self, theta, *args):
        """ Simple poisson log likelihood function
    
        :theta: model parameter (m500 [Msolar])
    
        """

        # standard slope parameter for beta model
        beta = 2./3
    
        # calc Lx from kettula15
        llx = self.calcLx(theta)[0]

        # calc r500 from mass estimate    
        lr500 = self.calcR5(theta)
    
        # fix beta model core radius (rc) to 1/3*r500
        rc = lr500 / 3.


        if self.aperture == None:
            # self-calibrate according to mass estimate in 0.2-0.5 r500 range
            regMask = np.logical_and(
                self.pixSep>=.2*lr500,
                self.pixSep<=.5*lr500
                )
        else:
            # use precomputed mask according to aperture values
            regMask = self.regMask


        # calculate observed counts in aperture
        obsCnts1D = self.imgData[regMask].flatten()
    
        # calculate Kettula+15 mask; Lx is calculated within 0.1-1 r500
        regMaskK15 = np.logical_and(
            self.pixSep>=.1*lr500,
            self.pixSep<=1.*lr500
            )
    
        # create model image, just takes the beta model value at the center of the pixel
        # it is not integrated over the pixel size
        # but that does not change much because it is smooth
        # implement oversampling and averaging?
        self.modelImg = self.singleBeta(self.pixSep, beta, rc)
    
        # normalize model to unity in K15 range
        self.modelImg /= np.sum(self.modelImg[regMaskK15].flatten())
        modelNorm = llx / self.lum_ctr
    
        self.modelImg *= modelNorm
        self.modelImg *= self.expData

        # add background counts to model image
        self.modelImg += self.bkgData
    
        # calculate model + background counts in aperture
        modelCnts1D = self.modelImg[regMask].flatten()
    
        # calculate total observed and model counts
        modelCnts = np.sum(modelCnts1D)
        obsCnts = np.sum(obsCnts1D)
        lnobsCntsFac = math.lgamma(obsCnts+1)

        return np.log(modelCnts)*obsCnts - modelCnts - lnobsCntsFac


    def lnprior(self, theta):    
        if (theta <= 0.): # mass must be positive
            return -np.inf
        else:
            return 0.


    def lnprob(self, theta, *args):
        lp = self.lnprior(theta)
        if not np.isfinite(lp):
          return -np.inf
        
        try:
          ll = self.lnlike(theta, *args)
        except Exception as excep:
          print(excep)
          return -np.inf

        if np.isnan(ll):
          return -np.inf

        return lp + ll


    def runFit(self, outDir, nwalkers=128, nstepsBurnin=500, nsteps=500):
        """ calculate mass
        """
        if not os.path.isdir(outDir):
            os.makedirs(outDir)

        # number of parameters; Is one, namely the cluster mass
        nparams = 1 
    
        # intial mass estimate
        m500t = 1e14 
    
        # randomize starting position for walkers
        p0 = [m500t]
        p0 = p0 + 0.1*m500t * np.random.randn(nwalkers, len(p0))

        # set up backend
        backend = emcee.backends.HDFBackend( os.path.join(outDir,'fit.h5') )
        backend.reset(nwalkers, nparams)

        samplerScatter = emcee.EnsembleSampler(nwalkers, nparams, self.lnprob, backend=backend)
        
        # run burn-in
        pos, _, _ = samplerScatter.run_mcmc(p0, nstepsBurnin, progress=True)
 
        samplerScatter.reset()

        # run mcmc, start at positions after burn-in
        for sample in samplerScatter.sample(pos, iterations=nsteps, progress=True):
            pass

        meanAccF = np.mean(samplerScatter.acceptance_fraction)
        print("Mean acceptance fraction: {0:.3f}".format(meanAccF))
   
        # calculate best-fit value
        samplesScatter = samplerScatter.get_chain(flat=True)
        self.m500 = self.calcPercentiles(*np.percentile(samplesScatter, [15.87, 50.,84.13]))

        # store best-fit value in csv file
        mCsv = os.path.join(outDir,'m500.csv')
        with open(mCsv, 'w') as file_:
            outStr = '{0:e},-{1:e},+{2:e}'.format(*self.m500)
            file_.write(outStr)

        print('M500: ({0:.2e} -{1:.2e} +{2:.2e}) * 10^14 Msolar'.format(*self.m500/1e14))