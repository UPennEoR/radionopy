'''
radiono.rm

purpose | Module used to create RM object

Classes
-------
RM | allows user to grab data from IONEX files by altaz or radec
'''
from __future__ import print_function
import os
import healpy as hp
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, Angle, Latitude, Longitude, AltAz, ICRS
import radiono as rad
from radiono import physics as phys, interp as itp, ionex_file as inx, utils
from itertools import groupby

class IonoMap(object):
    '''
    superobject to do altaz or radec RM map generation
    gets the RM over the sphere given IONEX input

    Methods
    -------
    radec | gets B_para, RM, and dRM arrays by using ra and dec
    _hp_arr | gets altitude and azimuth arrays for healpix map
    altaz | gets B_para, RM, and dRM by using alt and az
    to_map | writes data to map
    to_alm | writes data to alm
    '''
    def __init__(self, lat_str, lon_str, times, height=0, ionex_dir=rad.ionex_dir, rm_dir=rad.rm_dir):
        '''
        initalizes RM object

        Parameters
        ----------
        lat_str | str: latitude that the observation was taken at
        lon_str | str: longitude that the observation was taken at
        times | list[str]: list of dates and times for observation in form "YYYY-MM-DDTHH:MM:SS". The time is a sexagesimal number.
        height | Optional[float]: height observation taken at in meters
        nside | Optional[int]: resolution of rm map --defaults to 16
        ionex_dir | Optional[str]: directory in which ionex files are / will be located
        rm_dir | Optional[str]: directory in which RM data files are / will be located
        '''
        # if not isinstance(date_strs,(list,tuple)):
        #     raise TypeError('date_strs must be a list')

        self.lat_str = lat_str
        self.lon_str = lon_str
        self.times = times
        self.height = height
        self.ionex_dir = ionex_dir
        self.rm_dir = rm_dir
        self.nside = 16
        self.B_para = None
        self.RMs = None
        self.dRMs = None
        self.coordinates_flag = None

        # The input list of dates and times need not be sorted in order. So we
        # will sort them for internal use and return a dictionary of RM maps
        # that are keyed on the input list of input dates.

        # There may be an arbitrary number of days, each with an arbitrary set
        # of times on that day. For each day, all of the maps for the times on
        # that day will be computed from the same IONEX file, and we only want
        # to process each file once. Thus we must group the input dates by day,
        # and process all of the times in each group together.

        unique_days = [item.split(' ')[0] for item in self.times]
        self.unique_days = sorted(list(set(unique_days)))

        # a dictionary of lists of decimal times, keyed by each unique day
        # string in self.times
        self.day_groups = {}

        # note that the times strings are ISO formatted so sorting
        # lexiographically, which is what sorted() does, sorts them into
        # chronological order.
        test_function = lambda x: x.split(' ')[0]
        for key, group in groupby(sorted(self.times), test_function):
            self.day_groups[key] = list(group)

    @property
    def lat(self):
        if self.lat_str[-1]=='s': m=-1.
        else: m=1.
        return m*Latitude(Angle(self.lat_str[:-1]))

    @property
    def lon(self):
        if self.lon_str[-1]=='w': m=-1.
        else: m=1.
        return m*Longitude(Angle(self.lon_str[:-1]))

    @property
    def location(self):
        return EarthLocation(lat=self.lat, lon=self.lon, height=self.height * u.m)

    @property
    def npix(self):
        return hp.nside2npix(self.nside)

    def make_rm_dir(self, time_str,verbose=False):
        '''
        creates directory to hold rm files for a particular date

        Parameters
        ----------
        time_str | str: date of rms

        Returns
        -------
        str: rm directory path
        '''
        RM_dir = os.path.join(self.rm_dir, '{date}'.format(date=time_str))
        if not os.path.exists(RM_dir):
            if verbose: print('Creating directory %s'%RM_dir)
            os.makedirs(RM_dir)
        return RM_dir

    def ionex_data(self, year, month, day, ionex_dir=rad.ionex_dir, **kwargs):
        '''
        gathers all relevant IONEX info from file for specific date

        Parameters
        ----------
        year | int: year
        month | int: numbered month of the year
        day | int: numbered day of the month
        ionex_dir | Optional[str]: directory in which ionex files are / will be located
        verbose | Optional[bool]: whether to print values or not

        Returns
        -------
        tuple:
            array: tec healpix map
            array: rms tec healpix map
            float: ionosphere height in meters
        '''
        IONEX_file = inx.pull_IONEX_file(year, month, day)
        TEC, _, all_info = inx.get_IONEX_data(IONEX_file, **kwargs)

        tec_a, rms_a, ion_height = all_info[7:]

        # tec_hp = itp.ionex2healpix(tec_a, UT, TEC['lat'], TEC['lon'], **kwargs)
        # rms_hp = itp.ionex2healpix(rms_a, UT, TEC['lat'], TEC['lon'], **kwargs)
        #
        # return tec_hp, rms_hp, ion_height

        return tec_a, rms_a, ion_height, TEC

    def calc_radec_rm(self, ras, decs, verbose=False):
        #TODO: allow ra,dec to be floats, rather than arrays of floats
        # (to maintain single-pointing functionality)

        if not all((i<=2. * np.pi and i>=0.) for i in ras):
            raise ValueError('All RAs must be between 0 and 2*pi radians')
        if not all((i<=np.pi/2. and i>=-np.pi/2.) for i in decs):
            raise ValueError('All Decs must be between -pi/2 and pi/2 radians')

        self.coordinates_flag = 'J2000_RaDec'

        #final storage arrays
        # b_para_s = []
        # rm_s = []
        # drm_s = []
        # lsts_s = []
        # alt_src_s = []

        self.B_paras = {}
        self.RMs = {}
        self.dRMs = {}
        self.lsts_s = {}
        alt_src_s = {}

        for uday in self.day_groups:
            group = self.day_groups[uday]
            time_strs = [item.split(' ')[1] for item in group]

            UTs_dec = []
            for time_str in time_strs:
                Hr, Min, Sec = [float(x) for x in time_str.split(':')]
                dec_hour = Hr + Min / 60. + Sec / 3600.
                UTs_dec.append(dec_hour)

            year, month, day = [int(x) for x in uday.split('-')]

            #data aquisition
            tec_a, rms_a, ion_height, TEC = self.ionex_data(year, month, day)

            tec_hp = itp.ionex2healpix(tec_a, UTs_dec, TEC['lat'], TEC['lon'])
            rms_hp = itp.ionex2healpix(rms_a, UTs_dec, TEC['lat'], TEC['lon'])

            #temp storage arrays
            # alt_src_all = np.zeros([24,3072])
            lsts = []
            # RM_add = []
            # dRM_add = []

            RMs = []
            dRMs = []


            # predict the ionospheric RM for every hour within a day
            for ui,UT in enumerate(UTs_dec):
                time = Time(uday + ' ' +  time_strs[ui], format='iso', scale='utc') # XXX is this string reconstructed properly? not tested
                c_icrs = SkyCoord(ra=ras * u.radian, dec=decs * u.radian,
                                        location=self.location,
                                        obstime=time,
                                        frame='icrs')

                # Added to calculate LST for the given time
                c_local = AltAz(az=0.*u.deg,alt=90.*u.deg,obstime=time,location=self.location)
                c_local_Zeq = c_local.transform_to(ICRS)
                lsts.append(c_local_Zeq.ra.degree)

                # Transform given RA/Dec into alt/az
                c_altaz = c_icrs.transform_to('altaz')
                alt_src = np.array(c_altaz.alt.degree)
                az_src = np.array(c_altaz.az.degree)

                # AltAz.zen doesn't have method to return angle data
                zen_src = np.array(Angle(c_altaz.zen).degree)

                # Calculating the ion piercing point (IPP) depends on alt/az coords
                coord_lat, coord_lon, az_punct, zen_punct = phys.ipp(self.lat_str, self.lon_str, az_src, zen_src, ion_height)

                #XXX B_para calculated per UT
                #these are the data we care about
                B_para = phys.B_IGRF(year, month, day, coord_lat, coord_lon, ion_height, az_punct, zen_punct)
                TEC_path, RMS_TEC_path = itp.get_los_tec(tec_hp[ui], rms_hp[ui], coord_lat, coord_lon, zen_punct)
                RMs_ut = phys.RotationMeasure(TEC_path, B_para)
                dRMs_ut = phys.RotationMeasure(RMS_TEC_path, B_para)

                #TODO: replace append commands with numpy array indicies
                RMs.append(RMs_ut)
                dRMs.append(dRMs_ut)


            self.RMs.update({key:RMs[ti] for ti, key in enumerate(group)})
            self.dRMs.update({key:dRMs[ti] for ti, key in enumerate(group)})
            self.B_paras[uday] = B_para

    def calc_ionRIME_rm(self, verbose=False):

        hpxidx = np.arange(self.npix)
        theta, phi = hp.pix2ang(self.nside, hpxidx)
        R = hp.rotator.Rotator(rot=[0,-120.712])

        zen, az = R(theta, phi)

        # az = -az
        az[az < 0] += 2. * np.pi

        zen_src = np.degrees(zen)
        az_src = np.degrees(az)

        self.RMs = {}
        self.dRMs = {}
        self.B_paras = {}

        for uday in self.day_groups:

            group = self.day_groups[uday]
            time_strs = [item.split(' ')[1] for item in group]

            UTs_dec = []
            for time_str in time_strs:
                Hr, Min, Sec = [float(x) for x in time_str.split(':')]
                dec_hour = Hr + Min / 60. + Sec / 3600.
                UTs_dec.append(dec_hour)

            year, month, day = [int(x) for x in uday.split('-')]

            tec_a, rms_a, ion_height, TEC = self.ionex_data(year, month, day)

            tec_hp = itp.ionex2healpix(tec_a, UTs_dec, TEC['lat'], TEC['lon'])
            rms_hp = itp.ionex2healpix(rms_a, UTs_dec, TEC['lat'], TEC['lon'])


            coord_lat, coord_lon, az_punct, zen_punct = phys.ipp(self.lat_str, self.lon_str,
                                                                 az_src, zen_src,
                                                                 ion_height)
            #XXX B_para calculated per DAY
            B_para = phys.B_IGRF(year, month, day,
                                 coord_lat, coord_lon,
                                 ion_height, az_punct, zen_punct)

            RMs = []
            dRMs = []

            for ui,UT in enumerate(UTs_dec):

                TEC_path, RMS_TEC_path = itp.get_los_tec(tec_hp[ui], rms_hp[ui],
                                                          coord_lat, coord_lon,
                                                          zen_punct)

                RM_ut = phys.RotationMeasure(TEC_path, B_para)
                dRM_ut = phys.RotationMeasure(RMS_TEC_path, B_para)

                RMs.append(RM_ut)
                dRMs.append(dRM_ut)

            self.RMs.update({key:RMs[ti] for ti, key in enumerate(group)})
            self.dRMs.update({key:dRMs[ti] for ti, key in enumerate(group)})
            self.B_paras[uday] = B_para


    def make_radec_RM_maps(self):
        """
        Generates the full maps in ra/dec coordinates.
        """
        ra, dec = self._radec_arr()
        self.calc_radec_rm(ra, dec)

    def _radec_arr(self):
        '''
        generates array of ra's and dec's for healpix map

        Returns
        -------
        tuple:
            array[float]: array of ra coordinates in radians
            array[float]: array of dec coordinates in radians
        '''
        hpxidx = np.arange(self.npix)
        cza, phi = hp.pix2ang(self.nside, hpxidx)
        dec = np.pi/2. - cza
        phi_m = np.amax(phi) # because of discrete sampling this is  2*\pi - \eppsilon
        ra = phi_m -phi

        return ra, dec

    def _hp_arr(self):
        '''
        generates array of altitudes and azimuths for healpix map

        Returns
        -------
        tuple:
            array[float]: array of altitudes in degrees
            array[float]: array of azimuths in degrees
        '''
        ipix = np.arange(self.npix)
        theta, phi = hp.pix2ang(self.nside, ipix)

        alt_src = np.amax(theta) - theta
        az_src = phi

        return alt_src, az_src

    def altaz(self):
        '''
        outputs RM data from altaz calculation

        Returns
        -------
        tuple:
            array[float]: parallel B field array
            array[float]: RM data
            array[float]: RM error data
        '''
        # rm_s = []
        # drm_s = []
        # b_para_s = []

        self.RMs = {}
        self.dRMs = {}
        self.B_paras = {}

        alt_src, az_src = self._hp_arr()
        zen_src , _ = hp.pix2ang(self.nside, np.arange(self.npix))

        az_src = np.degrees(az_src)
        zen_src = np.degrees(zen_src)

        for uday in self.day_groups:

            group = self.day_groups[uday]
            time_strs = [item.split(' ')[1] for item in group]

            UTs_dec = []
            for time_str in time_strs:
                Hr, Min, Sec = [float(x) for x in time_str.split(':')]
                dec_hour = Hr + Min / 60. + Sec / 3600.
                UTs_dec.append(dec_hour)

            year, month, day = [int(x) for x in uday.split('-')]

            tec_a, rms_a, ion_height, TEC = self.ionex_data(year, month, day)

            tec_hp = itp.ionex2healpix(tec_a, UTs_dec, TEC['lat'], TEC['lon'])
            rms_hp = itp.ionex2healpix(rms_a, UTs_dec, TEC['lat'], TEC['lon'])


            coord_lat, coord_lon, az_punct, zen_punct = phys.ipp(self.lat_str, self.lon_str,
                                                                 az_src, zen_src,
                                                                 ion_height)
            #XXX B_para calculated per DAY
            B_para = phys.B_IGRF(year, month, day,
                                 coord_lat, coord_lon,
                                 ion_height, az_punct, zen_punct)

            RMs = []
            dRMs = []

            for ui,UT in enumerate(UTs_dec):

                TEC_path, RMS_TEC_path = itp.get_los_tec(tec_hp[ui], rms_hp[ui],
                                                          coord_lat, coord_lon,
                                                          zen_punct)

                RM_ut = phys.RotationMeasure(TEC_path, B_para)
                dRM_ut = phys.RotationMeasure(RMS_TEC_path, B_para)

                RMs.append(RM_ut)
                dRMs.append(dRM_ut)

            self.RMs.update({key:RMs[ti] for ti, key in enumerate(group)})
            self.dRMs.update({key:dRMs[ti] for ti, key in enumerate(group)})
            self.B_paras[uday] = B_para

def HERA_RM(times):
    """
    For our convenience: built-in generator for the PAPER/HERA site in the Karoo RQZ, South Africa
    """
    lat_str = '30d43m17.5ss' # -30.7215 degrees
    lon_str = '21d25m41.9se' # 21.4283 degrees
    height = 1073 # XXX
    return IonoMap(lat_str=lat_str, lon_str=lon_str, times=times, height=height)
