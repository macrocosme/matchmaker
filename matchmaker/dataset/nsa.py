import numpy as np
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
# from ..utils import load_fits_as_dataframe

h = 0.68 # Following Reines et al (2013) - used for mass

class Nsa(Catalog):
    file_location = DATA_BASE_PATH + 'data/NSA/nsa_v0_1_2.fits'
    name = 'nsa'

    def __init__(self, load_data=False, constrain=True, use_a=True, use_distance=False):
        super().__init__(ra=Column('RA', u.deg), dec=Column('DEC', u.deg), use_distance=use_distance)

        self.set_survey_specific_columns()

        self.cols.mstar = self.cols.all.mass
        self.cols.distance = self.cols.all.zdist
        self.cols.sfr = self.cols.all.b300 # Not sure which one to use yet
        self.cols.z = self.cols.all.z

        if load_data:
            self.load_data(constrain=constrain, use_a=use_a)

    @property
    def n_source(self):
        return len(self.df)

    def load_data(self, constrain=False, use_a=True):
        # This one requiires some specific trimming for our purpose
        from astropy.io import fits
        from astropy.table import Table

        with fits.open(self.file_location) as hdul:
            t = Table(hdul[1].data)

        t['ABSMAG F'] = t['ABSMAG'][:,0]
        t['ABSMAG N'] = t['ABSMAG'][:,1]
        t['ABSMAG u'] = t['ABSMAG'][:,2]
        t['ABSMAG g'] = t['ABSMAG'][:,3]
        t['ABSMAG r'] = t['ABSMAG'][:,4]
        t['ABSMAG i'] = t['ABSMAG'][:,5]
        t['ABSMAG z'] = t['ABSMAG'][:,6]
        t.remove_columns(['NMGY', 'NMGY_IVAR', 'RNMGY', 'ABSMAG', 'AMIVAR', 'EXTINCTION', 'KCORRECT', 'KCOEFF', 'MTOL',
                          'NPROF', 'PROFMEAN', 'PROFMEAN_IVAR', 'QSTOKES', 'USTOKES', 'BASTOKES', 'PHISTOKES',
                          'PETROFLUX', 'PETROFLUX_IVAR', 'FIBERFLUX', 'FIBERFLUX_IVAR', 'SERSICFLUX', 'SERSICFLUX_IVAR',
                          'ASYMMETRY', 'CLUMPY', 'DFLAGS', 'PROFTHETA', 'COEFF']) # For now... I'll come back to these multidimensional
                                                                                  # columns at some point.
        self.df = t.to_pandas()
        if constrain:
            df = self.df.loc[
                (self.df[self.cols.mstar.label] > 1e7) &
                (self.df[self.cols.mstar.label] <= 3e9) &
                (self.df[self.cols.all.absmag_g.label] > -20) &
                (self.df[self.cols.all.absmag_r.label] > -20)
            ]
            self.df = df

    def set_survey_specific_columns(self):
        # Radio columns
        self.cols.all = Column()
        self.cols.all.iauname = Column('IAUNAME', None,'IAU-style designation based on RA/Dec')
        self.cols.all.subdir = Column('SUBDIR', None,'Subdirectory for images in the detect directory')
        self.cols.all.ra = Column('RA', u.deg,'Right Ascension of measured object center (J2000)')
        self.cols.all.dec = Column('DEC', u.deg,'Declination of measured object center (J2000)')
        self.cols.all.isdss = Column('ISDSS', None,'Zero-indexed row of SDSS source file (-1 if no match)')
        self.cols.all.ined = Column('INED', None,'Zero-indexed row of NED source file (-1 if no match)')
        self.cols.all.isixdf = Column('ISIXDF', None,'Zero-indexed row of 6dFGRS source file (-1 if no match)')
        self.cols.all.ialfalfa = Column('IALFALFA', None,'Zero-indexed row of ALFALFA source file (-1 if no match)')
        self.cols.all.izcat = Column('IZCAT', None,'Zero-indexed row of ZCAT source file (-1 if no match)')
        self.cols.all.itwodf = Column('ITWODF', None,'Zero-indexed row of 2dFGRS source file (-1 if no match)')
        self.cols.all.mag = Column('MAG', u.mag,"'Nominal B magnitude; not reliable, only used to set size of image for analysis'")
        self.cols.all.z = Column('Z', None,'Heliocentric redshift')
        self.cols.all.zsrc = Column('ZSRC', None,"'Source of redshift determination (alfalfa, ned, sdss, sixdf, twodf, or zcat)'")
        self.cols.all.size = Column('SIZE', u.deg,'Size of analyzed mosaics')
        self.cols.all.run = Column('RUN', None,"'SDSS drift scan run covering catalog position (RACAT, DECCAT)'")
        self.cols.all.camcol = Column('CAMCOL', None,"'SDSS camcol covering catalog position (RACAT, DECCAT)'")
        self.cols.all.field = Column('FIELD', None,"'SDSS field covering catalog position (RACAT, DECCAT)'")
        self.cols.all.rerun = Column('RERUN', None,'Photometric rerun of SDSS used to determine SDSS coverage')
        self.cols.all.xpos = Column('XPOS', u.pix,"'SDSS camera column of catalog position (RACAT, DECCAT)'")
        self.cols.all.ypos = Column('YPOS', u.pix,"'SDSS camera row of catalog position (RACAT, DECCAT)'")
        self.cols.all.zlg = Column('ZLG', None,'Local Group redshift (following Yahil et al. 1977)')
        self.cols.all.zdist = Column('ZDIST', None, 'Distance estimate using pecular velocity model of Willick et al. (1997); multiply by c/H0 for Mpc')
        self.cols.all.zdist_err = Column('ZDIST_ERR', None,'Error estimate in distance')
        self.cols.all.nsaid = Column('NSAID', None,'unique ID within NSA catalog')
        self.cols.all.nmgy7 = Column('NMGY','nanomaggies','Galactic-extinction corrected AB flux used for K-correction (from SERSICFLUX) in FNugriz')
        self.cols.all.nmgy_ivar7 = Column('NMGY_IVAR','nanomaggies-2','Inverse variances used for K-correction in FNugriz (includes error floors)')
        self.cols.all.ok = Column('OK', None,"'1 if K-correction was performed, 0 if not (should all be 1)'")
        self.cols.all.rnmgy7 = Column('RNMGY','nanomaggies','Reconstructed AB nanomaggies from K-correction fit')
        self.cols.all.absmag_f = Column('ABSMAG F', u.mag,"'Absolute magnitude estimates for F from K-corrections (Ωm=0.3, ΩΛ=0.7, h=1)'")
        self.cols.all.absmag_n = Column('ABSMAG N', u.mag,"'Absolute magnitude estimates for N from K-corrections (Ωm=0.3, ΩΛ=0.7, h=1)'")
        self.cols.all.absmag_u = Column('ABSMAG u', u.mag,"'Absolute magnitude estimates for u from K-corrections (Ωm=0.3, ΩΛ=0.7, h=1)'")
        self.cols.all.absmag_g = Column('ABSMAG g', u.mag,"'Absolute magnitude estimates for g from K-corrections (Ωm=0.3, ΩΛ=0.7, h=1)'")
        self.cols.all.absmag_r = Column('ABSMAG r', u.mag,"'Absolute magnitude estimates for r from K-corrections (Ωm=0.3, ΩΛ=0.7, h=1)'")
        self.cols.all.absmag_i = Column('ABSMAG i', u.mag,"'Absolute magnitude estimates for i from K-corrections (Ωm=0.3, ΩΛ=0.7, h=1)'")
        self.cols.all.absmag_z = Column('ABSMAG z', u.mag,"'Absolute magnitude estimates for z from K-corrections (Ωm=0.3, ΩΛ=0.7, h=1)'")
        self.cols.all.amivar7 = Column('AMIVAR', u.mag**-2,'Inverse variance on absolute magnitudes in FNugriz')
        self.cols.all.extinction7 = Column('EXTINCTION', u.mag,"'Galactic extinction in FNugriz from Schlegel, Finkbeiner & Davis (1997)'")
        self.cols.all.kcorrect7 = Column('KCORRECT', u.mag,'K-corrections used for absolute magnitudes in FNugriz')
        self.cols.all.kcoeff5 = Column('KCOEFF[5]', None,'Coefficients of templates in K-correction fit (see Blanton & Roweis 2007)')
        self.cols.all.mtol7 = Column('MTOL',None,'Mass-to-light ratio derived in K-correction for FNugriz')
        self.cols.all.b300 = Column('B300', None,"'Star-formation rate b-parameter (current over past average) for last 300 Myrs (from K-correction fit, use with caution)'")
        self.cols.all.b1000 = Column('B1000', None,"'Star-formation rate b-parameter (current over past average) for last 1 Gyrs (from K-correction fit, use with caution)'")
        self.cols.all.mets = Column('METS', None,'Metallicity from K-correction fit (use with caution)')
        self.cols.all.mass = Column('MASS', u.Msun * h**-2,'Stellar mass from K-correction fit (use with caution)')
        self.cols.all.xcen = Column('XCEN', u.pix,'X centroid in child atlas image')
        self.cols.all.ycen = Column('YCEN', u.pix,'Y centroid in child atlas image')
        self.cols.all.nprof7 = Column('NPROF', None,'Number of measured azimuthally averaged profile bins in FNugriz')
        self.cols.all.profmean7_15 = Column('PROFMEAN','nanomaggies per pixel',"'Mean surface brightness in annuli in FNugriz (note that F and N have pixel size 1.5 u.arcsec, ugriz have pixel size 0.396 u.arcsec)'")
        self.cols.all.profmean_ivar7_15 = Column("'PROFMEAN_IVAR'",'nanomaggies-2','Inverse variance of PROFMEAN in FNugriz')
        self.cols.all.qstokes7_15 = Column("QSTOKES", None,'Q Stokes parameter in each annulus in FNugriz')
        self.cols.all.ustokes7_15 = Column("USTOKES", None,'U Stokes parameter in each annulus in FNugriz')
        self.cols.all.bastokes7_15 = Column("BASTOKES", None,'Axis ratio b/a derived from Stokes paramaters in each annulus in FNugriz')
        self.cols.all.phistokes7_15 = Column("PHISTOKES", u.deg,'Angle (E of N) derived from Stokes paramaters in each annulus in FNugriz')
        self.cols.all.petroflux7 = Column('PETROFLUX','nanomaggies','Azimuthally-averaged SDSS-style Petrosian flux in FNugriz (GALEX-SDSS photometric systems)')
        self.cols.all.petroflux_ivar7 = Column('PETROFLUX_IVAR','nanomaggies-2','Inverse variance of PETROFLUX')
        self.cols.all.fiberflux7 = Column('FIBERFLUX','nanomaggies','Flux in 3-u.arcsec diameter aperture (not apodized) in FNugriz')
        self.cols.all.fiberflux_ivar7 = Column('FIBERFLUX_IVAR','nanomaggies-2','Inverse variance of FIBERFLUX')
        self.cols.all.ba50 = Column('BA50', None,'Axis ratio b/a from Stokes parameters at 50% light radius')
        self.cols.all.phi50 = Column('PHI50', u.deg,'Angle (E of N) from Stokes parameters at 50% light radius')
        self.cols.all.ba90 = Column('BA90', None,'Axis ratio b/a from Stokes parameters at 90% light radius')
        self.cols.all.phi90 = Column('PHI90', u.deg,'Angle (E of N) from Stokes parameters at 90% light radius')
        self.cols.all.sersicflux7 = Column('SERSICFLUX','nanomaggies','2D Sersic fit flux in FNugriz (GALEX-SDSS photometric systems)')
        self.cols.all.sersicflux_ivar7 = Column('SERSICFLUX_IVAR','nanomaggies-2','Inverse variance of SERSICFLUX')
        self.cols.all.sersic_n = Column('SERSIC_N', None,'2D Sersic index from fit')
        self.cols.all.sersic_ba = Column('SERSIC_BA', None,'Axis ratio b/a from 2D Sersic fit')
        self.cols.all.sersic_phi = Column('SERSIC_PHI', u.deg,'Angle (E of N) of major axis in 2D Sersic fit')
        self.cols.all.asymmetry7 = Column('ASYMMETRY', None,'Asymetry parameter in FNugriz')
        self.cols.all.clumpy7 = Column('CLUMPY', None,'Clumpy parameter in FNugriz')
        self.cols.all.dflags7 = Column('DFLAGS', None,'Flags from measurement in each band (FNugriz)')
        self.cols.all.aid = Column('AID', None,'Child id within parent')
        self.cols.all.pid = Column('PID', None,'Parent id within mosaic for this object')
        self.cols.all.dversion = Column('DVERSION', None,'Version of dimage used for measurement')
        self.cols.all.proftheta15 = Column('PROFTHETA[15]', u.arcsec,'Outer radius of each measured annulus')
        self.cols.all.petrotheta = Column('PETROTHETA', u.arcsec,'Petrosian radius (derived from r band)')
        self.cols.all.petroth50 = Column('PETROTH50', u.arcsec,'Petrosian 50% light radius (derived from r band)')
        self.cols.all.petroth90 = Column('PETROTH90', u.arcsec,'Petrosian 90% light radius (derived from r band)')
        self.cols.all.sersic_th50 = Column('SERSIC_TH50', u.arcsec,'50% light radius of 2D Sersic fit (along major axis)')
        self.cols.all.objno = Column('OBJNO', None,'Index number for Renbin Yan fits to spectra (PLATE*1000+FIBERID)*100L+(MJD%100)')
        self.cols.all.plate = Column('PLATE', None,'SDSS plate number (0 if not observed)')
        self.cols.all.fiberid = Column('FIBERID', None,'SDSS fiber number (0 if not observed)')
        self.cols.all.mjd = Column('MJD', None,'SDSS MJD of spectroscopic observation (0 if not observed)')
        self.cols.all.coeff7 = Column('COEFF', None,'Coefficients for spectroscopic continuum fit')
        self.cols.all.vdisp = Column('VDISP', u.km*u.s**-1,'Velocity dispersion')
        self.cols.all.d4000 = Column('D4000', None,'Dn4000 parameter (Balogh et al 1999)')
        self.cols.all.d4000err = Column('D4000ERR', None,'Error in Dn4000 parameter')
        self.cols.all.fa = Column('FA', None,'A fraction in K+A fit')
        self.cols.all.faerr = Column('FAERR', None,'error in A fraction in K+A fit')
        self.cols.all.s2flux = Column('S2FLUX',1e-17 * u.erg * u.s**-1 * u.cm**-2,"'[SII] 6716, 6731 flux'")
        self.cols.all.s2fluxerr = Column('S2FLUXERR',1e-17 * u.erg * u.s**-1 * u.cm**-2,"'Error in [SII] 6716, 6731 flux'")
        self.cols.all.s2ew = Column('S2EW', None,"'[SII] 6716, 6731 equivalent width'")
        self.cols.all.s2ewerr = Column('S2EWERR', None,"'Error in [SII] 6716, 6731 equivalent width'")
        self.cols.all.s2vmeas = Column('S2VMEAS', u.km*u.s**-1,"'Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.s2vmerr = Column('S2VMERR', u.km*u.s**-1,"'Error in Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.s2ratio = Column('S2RATIO', None,'flux ratio in [SII] doublet ([SII] 6731 / [SII] 6716)')
        self.cols.all.haflux = Column('HAFLUX',1e-17 * u.erg * u.s**-1 * u.cm**-2,'Hα flux')
        self.cols.all.hafluxerr = Column('HAFLUXERR',1e-17 * u.erg * u.s**-1 * u.cm**-2,'Error in Hα flux')
        self.cols.all.haew = Column('HAEW', None,'Hα equivalent width')
        self.cols.all.haewerr = Column('HAEWERR', None,'Error in Hα equivalent width')
        self.cols.all.havmeas = Column('HAVMEAS', u.km*u.s**-1,"'Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.havmerr = Column('HAVMERR', u.km*u.s**-1,"'Error in Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.n2flux = Column('N2FLUX',1e-17 * u.erg * u.s**-1 * u.cm**-2,'[NII] 6584 flux')
        self.cols.all.n2fluxerr = Column('N2FLUXERR',1e-17 * u.erg * u.s**-1 * u.cm**-2,'Error in [NII] 6584 flux')
        self.cols.all.n2ew = Column('N2EW', None,'[NII] 6584 equivalent width')
        self.cols.all.n2ewerr = Column('N2EWERR', None,'Error in [NII] 6584 equivalent width')
        self.cols.all.n2vmeas = Column('N2VMEAS', u.km*u.s**-1,"'Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.n2vmerr = Column('N2VMERR', u.km*u.s**-1,"'Error in Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.hbflux = Column('HBFLUX',1e-17 * u.erg * u.s**-1 * u.cm**-2,'Hβ flux')
        self.cols.all.hbfluxerr = Column('HBFLUXERR',1e-17 * u.erg * u.s**-1 * u.cm**-2,'Error in Hβ flux')
        self.cols.all.hbew = Column('HBEW', None,'Hβ equivalent width')
        self.cols.all.hbewerr = Column('HBEWERR', None,'Error in Hβ equivalent width')
        self.cols.all.hbvmeas = Column('HBVMEAS', u.km*u.s**-1,"'Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.hbvmerr = Column('HBVMERR', u.km*u.s**-1,"'Error in Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.o1flux = Column('O1FLUX',1e-17 * u.erg * u.s**-1 * u.cm**-2,'[OI] 6300 flux')
        self.cols.all.o1fluxerr = Column('O1FLUXERR',1e-17 * u.erg * u.s**-1 * u.cm**-2,'Error in [OI] 6300 flux')
        self.cols.all.o1ew = Column('O1EW', None,'[OI] 6300 equivalent width')
        self.cols.all.o1ewerr = Column('O1EWERR', None,'Error in [OI] 6300 equivalent width')
        self.cols.all.o1vmeas = Column('O1VMEAS', u.km*u.s**-1,"'Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.o1vmerr = Column('O1VMERR', u.km*u.s**-1,"'Error in Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.o2flux = Column('O2FLUX',1e-17 * u.erg * u.s**-1 * u.cm**-2,'[OII] 3727 flux')
        self.cols.all.o2fluxerr = Column('O2FLUXERR',1e-17 * u.erg * u.s**-1 * u.cm**-2,'Error in [OII] 3727 flux')
        self.cols.all.o2ew = Column('O2EW', None,'[OII] 3727 equivalent width')
        self.cols.all.o2ewerr = Column('O2EWERR', None,'Error in [OII] 3727 equivalent width')
        self.cols.all.o2vmeas = Column('O2VMEAS', u.km*u.s**-1,"'Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.o2vmerr = Column('O2VMERR', u.km*u.s**-1,"'Error in Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.o3flux = Column('O3FLUX',1e-17 * u.erg * u.s**-1 * u.cm**-2,'[OIII] 5007 flux')
        self.cols.all.o3fluxerr = Column('O3FLUXERR',1e-17 * u.erg * u.s**-1 * u.cm**-2,'Error in [OIII] 5007 flux')
        self.cols.all.o3ew = Column('O3EW', None,'[OIII] 5007 equivalent width')
        self.cols.all.o3ewerr = Column('O3EWERR', None,'Error in [OIII] 5007 equivalent width')
        self.cols.all.o3vmeas = Column('O3VMEAS', u.km*u.s**-1,"'Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.o3vmerr = Column('O3VMERR', u.km*u.s**-1,"'Error in Gaussian σ for this line, corrected for instrumental broadening'")
        self.cols.all.ahgew = Column('AHGEW', None,'Hγ absorption equivalent width (no emission infill correction)')
        self.cols.all.ahgewerr = Column('AHGEWERR', None,'Error in Hγ absorption equivalent width')
        self.cols.all.ahdew = Column('AHDEW', None,'Hδ absorption equivalent width (no emission infill correction)')
        self.cols.all.ahdewerr = Column('AHDEWERR', None,'Error in Hδ absorption equivalent width')
        self.cols.all.ne3ew = Column('NE3EW', None,'[NeIII] equivalent width')
        self.cols.all.ne3ewerr = Column('NE3EWERR', None,'Error in [NeIII] equivalent width')
        self.cols.all.ne5ew = Column('NE5EW', None,'[NeV] equivalent width')
        self.cols.all.ne5ewerr = Column('NE5EWERR', None,'Error in [NeV] equivalent width')
        self.cols.all.av = Column('AV', None,"'AV assuming intrinsic Balmer decrement of 2.85, RV = 3.1 and O'Donnell (1994) extinction'")
        self.cols.all.s2nsamp = Column('S2NSAMP', None,'median signal-to-noise ratio per pixel between 3600 and 6900 A (restframe)')
        self.cols.all.racat = Column('RACAT', u.deg,'Right Ascension of catalog object (J2000)')
        self.cols.all.deccat = Column('DECCAT', u.deg,'Declination of catalog object (J2000)')
        self.cols.all.zsdssline = Column('ZSDSSLINE', None,'Redshift from emission lines')
        self.cols.all.survey = Column('SURVEY', None,'Survey within SDSS of plate')
        self.cols.all.programname = Column('PROGRAMNAME', None,'Program name of plate')
        self.cols.all.platequality = Column('PLATEQUALITY', None,'Quality of plate')
        self.cols.all.tile = Column('TILE', None,'Tile number of spectrum')
        self.cols.all.plug_ra = Column('PLUG_RA', u.deg,'Right Ascension of spectroscopic fiber (J2000)')
        self.cols.all.plug_dec = Column('PLUG_DEC', u.deg,'Declination of spectroscopic fiber (J2000)')
