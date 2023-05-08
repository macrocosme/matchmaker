import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)


class Vlass(Catalog):
    file_location = DATA_BASE_PATH + 'vlass/vlass_table1_v1_good_partial.csv'
    name = 'vlass'


    def __init__(self, load_data=False, constrain=True):
        super().__init__(ra=Column('RA', u.deg), dec=Column('DEC', u.deg))
        self.cols.ra_err = Column('E_RA', u.deg)
        self.cols.dec_err = Column('E_DEC', u.deg)
        self.cols.major = Column('DC_Maj', u.arcsec)
        self.cols.minor = Column('DC_Min', u.arcsec)
        self.cols.pa = Column('DC_PA', u.deg)

        self.beam = 2.5 * u.arcsec
        self.frequency_band = 'S'
        self.central_frequency = 3000 * u.MHz

        self.meta = {'observation_period': '2017::2019'}

        self.set_survey_specific_columns()
        self.cols.measure = Column()
        self.cols.measure.total_flux = self.cols.component.flux.total_flux
        self.cols.measure.total_flux_err = self.cols.component.flux.e_total_flux
        self.cols.measure.peak_flux = self.cols.component.flux.peak_flux
        self.cols.measure.peak_flux_err = self.cols.component.flux.e_peak_flux

        if load_data:
            self.load_data()

    def load_data(self):
        df = pd.read_csv(self.file_location)

        # Some RA and DEC uncertainties are marked as -99
        # which makes the code crash later, so I replace it with 0 for now.
        # The following comparison may break at some point...
        df.loc[df[self.cols.ra_err.label] == -99] = 0.
        df.loc[df[self.cols.dec_err.label] == -99] = 0.

        self.df = df

    # Should simply fetch the same code from the Lotss class as it is identical.
    def semi_major(self, mask=None, with_unit=False, to_unit=None):
        a = self.prop_to_unit('major', self.cols.major.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                a = (a * self.cols.major.unit).to(to_unit).value
            else:
                a = a.to(to_unit)

        return a

    # idem
    def semi_minor(self, mask=None, with_unit=False, to_unit=None):
        b = self.prop_to_unit('minor', self.cols.minor.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                b = (b * self.cols.minor.unit).to(to_unit).value
            else:
                b = b.to(to_unit)

        return b

    def set_survey_specific_columns(self):
        self.cols.component = Column()
        self.cols.component.component_name = Column('Component_name', "Unique name of the PyBDSF component")
        self.cols.component.component_id = Column('Component_id', "The ID of the PyBDSF assigned component (referred to as ‘sources’ by PyBDSF) - unique only within subtile.")

        self.cols.component.coords = Column()
        self.cols.component.coords.ra = Column('RA', u.deg, "Right ascension of component")
        self.cols.component.coords.dec = Column('DEC', u.deg, "Declination of component")
        self.cols.component.coords.e_ra = Column('E_RA', u.deg, "Error in RA")
        self.cols.component.coords.e_dec = Column('E_DEC', u.deg, "Error in dec")

        self.cols.component.flux = Column()
        self.cols.component.flux.total_flux = Column('Total_flux', u.mJy, "Integrated flux of the component")
        self.cols.component.flux.e_total_flux = Column('E_Total_flux', u.mJy, "Error in integrated flux")
        self.cols.component.flux.peak_flux = Column('Peak_flux', u.mJy * u.beam**-1, "Peak flux of component")
        self.cols.component.flux.e_peak_flux = Column('E_Peak_flux', u.mJy * u.beam**-1, "Error in component peak flux")

        self.cols.component.isl_id = Column('Isl_id', "The ID of the PyBDSF island in which the component is located - unique only within subtile.")
        self.cols.component.maj = Column('Maj', u.arcsec, "Component major axis size (FWHM)")
        self.cols.component.e_maj = Column('E_Maj', u.arcsec, "Error in Maj")
        self.cols.component.min = Column('Min', u.arcsec, "Component minor axis size (FWHM)")
        self.cols.component.e_min = Column('E_Min', u.arcsec, "Error in Min")
        self.cols.component.pa = Column('PA', u.deg, "Component position angle, east of north")
        self.cols.component.e_pa = Column('E_PA', u.deg, "Error in PA")
        self.cols.component.isl_total_flux = Column('Isl_Total_flux', u.mJy, "Total flux in island containing the component")
        self.cols.component.e_isl_total_flux = Column('E_Isl_Total_flux', u.mJy, "Error in Isl total flux")
        self.cols.component.isl_rms = Column('Isl_rms', u.mJy * u.beam**-1, "Local rms of island containing component")
        self.cols.component.isl_mean = Column('Isl_mean', u.mJy * u.beam**-1, "Mean background flux of island containing the component")
        self.cols.component.resid_isl_rms = Column('Resid_Isl_rms', u.mJy * u.beam**-1, "Average residual background rms of the island")
        self.cols.component.resid_isl_mean = Column('Resid_Isl_mean', u.mJy * u.beam**-1, "Average residual background mean flux of Island")
        self.cols.component.ra_max = Column('RA_max', u.deg, "RA of the maximum flux of the component (peak flux position")
        self.cols.component.dec_max = Column('DEC_max', u.deg, "Dec of the maximum flux of the component (peak flux position)")
        self.cols.component.e_ra_max = Column('E_RA max', u.deg, "1 error in RA max")
        self.cols.component.e_dec_max = Column('E_DEC max', u.deg, "1 error in DEC max")
        self.cols.component.s_code = Column('S_code', None, "PyBDSF component type")
        self.cols.component.xposn = Column('Xposn', u.pixel, "X-position of component in the image array")
        self.cols.component.e_xposn = Column('E_Xposn', u.pixel, "1 error in Xposn")
        self.cols.component.yposn = Column('Yposn', u.pixel, "Y-position of component in the image array")
        self.cols.component.e_yposn = Column('E_Yposn', u.pixel, "1 error in Yposn")
        self.cols.component.xposn_max = Column('Xposn_max', u.pixel, "X-position of component peak in the image array")
        self.cols.component.e_xposn_max = Column('E_Xposn_max', u.pixel, "1 error in Xposn max")
        self.cols.component.yposn_max = Column('Yposn_max', u.pixel, "Y-position of component peak in the image array")
        self.cols.component.e_yposn_max = Column('E_Yposn_max', u.pixel, "1error in Yposn max")
        self.cols.component.maj_img_plane = Column('Maj_img_plane', u.arcsec, "Component major axis size in the image plane (FWHM)")
        self.cols.component.e_maj_img_plane = Column('E_Maj_img_plane', u.arcsec, "1error in E Maj img plane")
        self.cols.component.min_img_plane = Column('Min_img_plane', u.arcsec, "Component minor axis size in the image plane (FWHM)")
        self.cols.component.e_min_img_plane = Column('E_Min_img_plane', u.arcsec, "1error in Min img plane")
        self.cols.component.pa_img_plane = Column('PA_img_plane', u.deg, "Component position angle in the image plane, east of north")
        self.cols.component.e_pa_img_plane = Column('E_PA img plane', u.deg, "1error in PA img plane")
        self.cols.component.dc_maj = Column('DC_Maj', u.arcsec, "Deconvolved component major axis size(FWHM)")
        self.cols.component.e_dc_maj = Column('E_DC_Maj', u.arcsec, "1error in DC Maj")
        self.cols.component.dc_min = Column('DC_Min', u.arcsec, "Deconvolved component minor axis size(FWHM)")
        self.cols.component.e_dc_min = Column('E_DC_Min', u.arcsec, "1error in DC Min")
        self.cols.component.dc_pa = Column('DC_PA', u.deg, "Deconvolved component position angle, east of north")
        self.cols.component.e_dc_pa = Column('E_DC_PA', u.deg, "1error in DC PA")
        self.cols.component.dc_maj_img_plane = Column('DC_Maj_img_plane', u.arcsec, "Deconvolved component major axis size in the image plane (FWHM)")
        self.cols.component.e_dc_maj_img_plane = Column('E_DC_Maj_img_plane', u.arcsec, "1error in DC Maj img plane")
        self.cols.component.dc_min_img_plane = Column('DC_Min_img_plane', u.arcsec, "Deconvolved component minor axis size in the image plane (FWHM)")
        self.cols.component.e_dc_min_img_plane = Column('E_DC_Min_img_plane', u.arcsec, "1error in DC Min img plane")
        self.cols.component.dc_pa_img_plane = Column('DC_PA_img_plane', u.deg, "Deconvolved component position angle in the image plane, east of north")
        self.cols.component.e_dc_pa_img_plane = Column('E_DC_PA_img_plane', u.deg, "1error in DC PA img plane")
        self.cols.component.tile = Column('Tile', None, "VLASS tile this component is located in")
        self.cols.component.subtile = Column('Subtile', None, "VLASS subtile this component is located in")

        self.cols.distance = Column()
        self.cols.distance.nvss_distance = Column('NVSS_distance', u.arcsec, "Distance to nearest NVSS component")
        self.cols.distance.first_distance = Column('FIRST_distance', u.arcsec, "Distance to nearest FIRST component.")

        self.cols.other = Column()
        self.cols.other.peak_to_ring = Column('Peak_to_ring', None, "Ratio of the Peak flux to the maximum flux in annulus of r, centred on component RA, Dec.")
        self.cols.other.duplicate_flag = Column('Duplicate_flag', None, "Flag to denote the duplicate status of a component")
        self.cols.other.quality_flag = Column('quality_flag', None, "QA flag to deal with spurious detections and duplicates due to overlap between VLASS tile edges.")
        self.cols.other.ql_image_ra = Column('QL_image_RA', u.deg, "Central RA of the the QL image the source is in")
        self.cols.other.ql_image_dec = Column('QL_image_Dec', u.deg, "Central Dec of the the QL image the source is in")
        self.cols.other.ql_cutout = Column('QL_cutout', None, "url to 3⇥30 cutout of VLASS QL image centred on component position.")

        self.cols.source = Column()
        self.cols.source.source_name = Column('Source_name', None, "Name of the distinct radio source to which this component has been assigned")
        self.cols.source.source_type = Column('Source_type', None, "Flag to indicate the type of source the component has been assigned to")
        self.cols.component.nn_dist = Column('NN_dist', u.arcsec, "Angular distance to nearest other component in the Component Table")

        self.cols.beam = Column()
        self.cols.beam.bmaj = Column('BMAJ', u.arcsec, "Major axis of the beam")
        self.cols.beam.bmin = Column('BMIN', u.arcsec, "Minor axis of the beam")
        self.cols.beam.bpa = Column('BPA', u.deg, "Beam position angle, east of north")
