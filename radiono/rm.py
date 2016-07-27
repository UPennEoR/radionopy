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
from astropy.coordinates import SkyCoord, EarthLocation, Angle, Latitude, Longitude
import radiono as rad
from radiono import physics as phys, interp as itp, ionex_file as inx

class RM(object):
    '''
    superobject to do altaz or radec RM value generation
    gets the RM based on input

    Methods
    -------
    radec | gets B_para, RM, and dRM arrays by using ra and dec
    _hp_arr | gets altitude and azimuth arrays for healpix map
    altaz | gets B_para, RM, and dRM by using alt and az
    to_map | writes data to map
    to_alm | writes data to alm
    '''
    def __init__(self, lat_str, lon_str, time_strs, height=0, nside=16, ionex_dir=rad.ionex_dir, rm_dir=rad.rm_dir):
        '''
        initalizes RM object

        Parameters
        ----------
        lat_str | str: latitude that the observation was taken at
        lon_str | str: longitude that the observation was taken at
        time_strs | list[str]: list of dates for observation
        height | Optional[float]: height observation taken at in meters
        nside | Optional[int]: resolution of rm map --defaults to 16
        ionex_dir | Optional[str]: directory in which ionex files are / will be located
        rm_dir | Optional[str]: directory in which RM data files are / will be located
        '''
        self.lat_str = lat_str
        self.lon_str = lon_str
        self.times = Time(time_strs, format='isot')
        self.height = height
        self.ionex_dir = ionex_dir
        self.rm_dir = rm_dir
        self.nside = nside
        self.B_para = None
        self.RMs = None
        self.dRMs = None
        self.UTs = np.linspace(0, 23, num=24)


    @property
    def lat(self):
        return Latitude(Angle(self.lat_str[:-1]))

    @property
    def lon(self):
        return Longitude(Angle(self.lon_str[:-1]))

    @property
    def location(self):
        return EarthLocation(lat=self.lat, lon=self.lon, height=self.height * u.m)

    def make_rm_dir(self, time_str):
        '''
        creates directory to hold rm files for a particular date

        Parameters
        ----------
        time_str | str: date of rms

        Returns
        -------
        str: rm directory path
        '''
        RM_dir = os.path.join(self.rm_dir, '{date}'.format(date=time_str.split('T')[0]))
        if not os.path.exists(RM_dir):
            os.makedirs(RM_dir)
        return RM_dir

    def ionex_data(self, year, month, day, ionex_dir=rad.ionex_dir, verbose=False):
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
        IONEX_file = inx.IONEX_file_needed(year, month, day)
        TEC, _, all_info = inx.read_IONEX_TEC(IONEX_file, verbose=verbose)

        tec_a, rms_a, ion_height = all_info[7:]

        tec_hp = itp.interp_time(tec_a, TEC['lat'], TEC['lon'], verbose=verbose)
        rms_hp = itp.interp_time(rms_a, TEC['lat'], TEC['lon'], verbose=verbose)

        ## an idea, to have interp_time give maps at an arbitrary number of times throughout the day.
        ## Not yet developed.
        # tec_hp = itp._interp_time(tec_a, TEC['lat'], TEC['lon'], ntimes=ntimes, verbose=verbose)
        # rms_hp = itp._interp_time(rms_a, TEC['lat'], TEC['lon'], ntimes=ntimes, verbose=verbose)

        return tec_hp, rms_hp, ion_height

    def radec(self, ras, decs):
        '''
        outputs RM data from radec calculation

        Parameters
        ----------
        ras | array[float]: array of RAs for observation, corresponds by element to decs
        decs | array[float]: array of DECs for observation, corresponds by element for ras
        '''
        for time in self.times:
            time_str = str(time)
            RM_dir = self.make_rm_dir(time_str)

            year, month, day = time_str.split('T')[0].split('-')
            tec_hp, rms_hp, ion_height = self.ionex_data(year, month, day)

            # predict the ionospheric RM for every hour within a day 
            for UT in self.UTs:
                hour = rad.std_hour(UT)
                ra_dec = SkyCoord(ra=ras, dec=decs,
                                  location=self.location, obstime=time + UT * u.hr)
                altaz = ra_dec.altaz

                alt_src = altaz.alt
                az_src = altaz.az
                zen_src = altaz.zen

                if len(alt_src.shape) <= 1:
                    #alt_src = np.array([alt_src.item().degree])
                    #az_src = np.array([az_src.item().degree])
                    #zen_src = np.array([zen_src.value])
                    alt_src = alt_src.item().degree
                    az_src = az_src.item().degree
                    zen_src = zen_src.value

                coord_lat, coord_lon,\
                az_punct, zen_punct = phys.ipp(self.lat_str, self.lon_str,
                                               az_src, zen_src, ion_height)

                B_para = phys.B_IGRF(year, month, day,
                                     coord_lat, coord_lon,
                                     ion_height, az_punct, zen_punct)

                TEC_path, RMS_TEC_path = itp.interp_space(tec_hp[UT], rms_hp[UT],
                                                          coord_lat, coord_lon,
                                                          zen_punct)

                new_file = os.path.join(RM_dir, 'IonRM{hour}.txt'.format(hour=hour))
                rad.write_RM(hour, new_file, B_para, TEC_path, RMS_TEC_path, write_to_file=True)
        self.parse_RM()

    def parse_RM(self):
        '''
        parses ionospheric RM files and assigns values to object

        Returns
        -------
        tuple:
            array[float]: parallel B field array
            array[float]: RM data
            array[float]: RM error data
        '''
        b_para_s = []
        rm_s = []
        drm_s = []
        for time in self.times:
            time_str = str(time)
            RM_add = []
            dRM_add = []
            RM_dir = os.path.join(self.rm_dir, '{date}'.format(date=time_str.split('T')[0]))
            for UT in self.UTs:
                data_file = os.path.join(RM_dir, 'IonRM{hour}.txt'.format(hour=rad.std_hour(UT, verbose=False)))
                _, _, B_para, RM_ut, dRM_ut = np.loadtxt(data_file, unpack=True)
                b_para_s.append(B_para)
                RM_add.append(RM_ut)
                dRM_add.append(dRM_ut)
            rm_s.append(RM_add)
            drm_s.append(dRM_add)

        self.B_para = np.array(b_para_s)
        self.RMs = np.array(rm_s)
        self.dRMs = np.array(drm_s)

        return self.B_para, self.RMs, self.dRMs

    def _hp_arr(self):
        '''
        generates array of altitudes and azimuths for healpix map

        Returns
        -------
        tuple:
            array[float]: array of altitudes
            array[float]: array of azimuths
        '''
        npix = hp.nside2npix(self.nside)
        ipix = np.arange(npix)
        theta, phi = hp.pix2ang(self.nside, ipix)

        alt_src = 90. - np.degrees(np.array(theta))
        az_src = np.degrees(np.array(phi))

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
        rm_s = []
        drm_s = []
        b_para_s = []

        alt_src, az_src = self._hp_arr()
        zen_src = 90. - alt_src

        for date in self.times:
            date_str = str(date)
            RM_dir = self.make_rm_dir(date_str)
            year, month, day = date_str.split('T')[0].split('-')
            tec_hp, rms_hp, ion_height = self.ionex_data(year, month, day)

            coord_lat, coord_lon, az_punct, zen_punct = phys.ipp(self.lat_str, self.lon_str,
                                                                 az_src, zen_src,
                                                                 ion_height)

            B_para = phys.B_IGRF(year, month, day,
                                 coord_lat, coord_lon,
                                 ion_height, az_punct, zen_punct)

            RMs = []
            dRMs = []
            for UT in self.UTs:
                hour = rad.std_hour(UT)
                radec_file = os.path.join(RM_dir, 'radec{hour}.txt'.format(hour=hour))
                rad.write_radec(UT, radec_file, alt_src, az_src, date_str, self.location)

                TEC_path, RMS_TEC_path = itp.interp_space(tec_hp[UT], rms_hp[UT],
                                                          coord_lat, coord_lon,
                                                          zen_punct)

                new_file = os.path.join(RM_dir, 'IonRM{hour}.txt'.format(hour=hour))
                rad.write_RM(hour, new_file, B_para, TEC_path, RMS_TEC_path)

                _, _, _, RM_add, dRM_add = np.loadtxt(new_file, unpack=True)
                RMs.append(RM_add)
                dRMs.append(dRM_add)

            rm_s.append(RM_add)
            drm_s.append(dRM_add)
            b_para_s.append(B_para)

        self.RMs = np.array(rm_s)
        self.dRMs = np.array(drm_s)
        self.B_paras = np.array(b_para_s)

    def to_map(self, filename):
        '''
        writes map from RM data

        Parameters
        ----------
        filename | str: path to file to write to
        '''
        if self.RMs is None:
            raise Exception
        hp.fitsfunc.write_map(filename, self.RMs)


    def to_alm(self, map_file, alm_file):
        '''
        writes alm from map file

        Parameters
        ----------
        map_file | str: map file path
        alm_file | str: path to file to write to
        '''
        alms = hp.sphtfunc.map2alm(map_file)
        hp.fitsfunc.write_alm(alm_file, alms)

    def map_to_npz(self, time_str, loc_str='PAPER', verbose=True):
        '''
        writes map of particular time to npz file

        Parameters
        ----------
        time_str | str: time and date of map
        loc_str | str: name of location to incorporate into filename output
        verbose | Optional[bool]: whether to print values or not
        '''
        #I could fish around in the file read to get npix and then re-loop, but why not just be lazy sometimes
        rng = np.arange(24)
        final_TEC, final_rm, final_drm, ra, dec = np.zeros((rng.shape[0], self.npix)),\
                                                  np.zeros((rng.shape[0], self.npix)),\
                                                  np.zeros((rng.shape[0], self.npix)),\
                                                  np.zeros((rng.shape[0], self.npix)),\
                                                  np.zeros((rng.shape[0], self.npix))
        RM_dir = os.path.join(self.rm_dir, '{date}'.format(date=time_str.split('T')[0]))
        for UT in rng:
            rm_file = os.path.join(RM_dir, 'IonRM{num}.txt'.format(num=rad.std_hour(UT, verbose=verbose)))
            radec_file = os.path.join(RM_dir, 'radec{num}.txt'.format(num=rad.std_hour(UT, verbose=verbose)))
            
            _, TEC, B, RM, dRM = np.loadtxt(rm_file, unpack=True)
            RA, DEC = np.loadtxt(radec_file, unpack=True)
            
            final_TEC[UT, :] = TEC
            final_rm[UT, :] = RM
            final_drm[UT, :] = dRM
            ra[UT,:] = RA
            dec[UT,:] = DEC

        f_name = ''.join((time_str.split('T')[0], '_', loc_str, '.npz'))
        npz_dir = os.path.join(rad.root_dir, 'npz')
        if not os.path.exists(npz_dir):
            os.mkdir(npz_dir)
        npz_file = os.path.join(npz_dir, f_name)
        if verbose:
            print('Saving TEC, RM +/- dRM data and RA/Dec mapping to {filename}'.format(filename=npz_file))

        np.savez(npz_file, TEC=final_TEC, RM=final_rm, dRM=final_drm, RA=ra, DEC=dec)

    def maps_to_npz(self):
        '''
        writes all maps to npz files
        '''
        for time in self.times:
            time_str = str(time)
            self.map_to_npz(time_str)

if __name__ == '__main__':
    print('This is not a script anymore')
