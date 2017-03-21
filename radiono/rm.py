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
    def __init__(self, lat_str, lon_str, date_strs, height=0, ionex_dir=rad.ionex_dir, rm_dir=rad.rm_dir):
        '''
        initalizes RM object

        Parameters
        ----------
        lat_str | str: latitude that the observation was taken at
        lon_str | str: longitude that the observation was taken at
        date_strs | list[str]: list of dates for observation in form "YYYY-MM-DD"
        height | Optional[float]: height observation taken at in meters
        nside | Optional[int]: resolution of rm map --defaults to 16
        ionex_dir | Optional[str]: directory in which ionex files are / will be located
        rm_dir | Optional[str]: directory in which RM data files are / will be located
        '''
        if not isinstance(date_strs,(list,tuple)):
            raise TypeError('date_strs must be a list')

        self.lat_str = lat_str
        self.lon_str = lon_str
        self.times = Time(date_strs, format='isot')
        self.height = height
        self.ionex_dir = ionex_dir
        self.rm_dir = rm_dir
        self.nside = 16
        self.B_para = None
        self.RMs = None
        self.dRMs = None
        self.UTs = np.linspace(0, 23, num=24)
        self.coordinates_flag = None


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

        tec_hp = itp.ionex2healpix(tec_a, TEC['lat'], TEC['lon'], **kwargs)
        rms_hp = itp.ionex2healpix(rms_a, TEC['lat'], TEC['lon'], **kwargs)

        return tec_hp, rms_hp, ion_height

    def calc_radec_rm(self, ras, decs, verbose=False):
        #TODO: allow ra,dec to be floats, rather than arrays of floats 
        # (to maintain single-pointing functionality)
         
        if not all((i<=2. * np.pi and i>=0.) for i in ras):
            raise ValueError('All RAs must be between 0 and 2*pi radians')
        if not all((i<=np.pi/2. and i>=-np.pi/2.) for i in decs):
            raise ValueError('All Decs must be between -pi/2 and pi/2 radians')

        self.coordinates_flag = 'J2000_RaDec'
        
        #final storage arrays
        b_para_s = []
        rm_s = []
        drm_s = []
        lsts_s = []
        alt_src_s = []
        
        for time in self.times: #for every day
            #parsing, creating directory structure
            time_str,_,_ = str(time).partition('T')
            RM_dir = self.make_rm_dir(time_str)
            year, month, day = map(int,time_str.split('-'))
            if verbose: print(year,month,day)
            
            #data aquisition
            tec_hp, rms_hp, ion_height = self.ionex_data(year, month, day)
           
            #temp storage arrays
            alt_src_all = np.zeros([24,3072])
            lsts = np.zeros(24)
            RM_add = []
            dRM_add = []            
            # predict the ionospheric RM for every hour within a day
            for ui,UT in enumerate(self.UTs):
                hour = utils.std_hour(UT)
                c_icrs = SkyCoord(ra=ras * u.radian, dec=decs * u.radian,
                                        location=self.location, obstime=time + UT * u.hr, frame='icrs')

                # Added to calculate LST for the given time
                c_local = AltAz(az=0.*u.deg,alt=90.*u.deg,obstime=time + UT * u.hr,location=self.location)
                c_local_Zeq = c_local.transform_to(ICRS)
                lsts[ui] = c_local_Zeq.ra.degree

                # Transform given RA/Dec into alt/az
                c_altaz = c_icrs.transform_to('altaz')
                alt_src = np.array(c_altaz.alt.degree)
                az_src = np.array(c_altaz.az.degree)
                alt_src_all[ui,:] = alt_src
                
                # AltAz.zen doesn't have method to return angle data
                zen_src = np.array(Angle(c_altaz.zen).degree)
                
                # Calculating the ion piercing point (IPP) depends on alt/az coords
                coord_lat, coord_lon, az_punct, zen_punct = phys.ipp(self.lat_str, self.lon_str, az_src, zen_src, ion_height)
                
                #XXX B_para calculated per UT
                #these are the data we care about
                B_para = phys.B_IGRF(year, month, day, coord_lat, coord_lon, ion_height, az_punct, zen_punct)
                TEC_path, RMS_TEC_path = itp.get_los_tec(tec_hp[ui], rms_hp[ui], coord_lat, coord_lon, zen_punct)
                IRM = 2.6e-17 * B_para * TEC_path
                rms_IRM = 2.6e-17 * B_para * RMS_TEC_path
                
                #TODO: replace append commands with numpy array indicies               
                b_para_s.append(B_para)
                RM_add.append(IRM)
                dRM_add.append(rms_IRM)
            
            alt_src_s.append(alt_src_all)
            lsts_s.append(lsts)    
            rm_s.append(RM_add)
            drm_s.append(dRM_add)
            
        self.B_para = np.array(b_para_s)
        self.RMs = np.array(rm_s)
        self.dRMs = np.array(drm_s)
        self.alt_src = np.array(alt_src_s)
        self.lst = np.array(lsts_s)

    def make_radec_RM_maps(self):
        """
        Generates the full maps in ra/dec coordinates.
        """
        ra, dec = self._radec_arr()
        self.calc_radec_RM(ra, dec)

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
            date_str,_,_ = str(date).partition('T')
            RM_dir = self.make_rm_dir(date_str)
            year, month, day = map(int,date_str.split('-'))
            tec_hp, rms_hp, ion_height = self.ionex_data(year, month, day)

            coord_lat, coord_lon, az_punct, zen_punct = phys.ipp(self.lat_str, self.lon_str,
                                                                 az_src, zen_src,
                                                                 ion_height)
            #XXX B_para calculated per DAY
            B_para = phys.B_IGRF(year, month, day,
                                 coord_lat, coord_lon,
                                 ion_height, az_punct, zen_punct)

            RMs = []
            dRMs = []
            for ui,UT in enumerate(self.UTs):
                hour = utils.std_hour(UT)
                radec_file = os.path.join(RM_dir, 'radec{hour}.txt'.format(hour=hour))
                utils.write_radec(UT, radec_file, alt_src, az_src, date_str, self.location)
                #UTs are integer distances apart. Would an enumeration be better?
                TEC_path, RMS_TEC_path = itp.get_los_tec(tec_hp[ui], rms_hp[ui],
                                                          coord_lat, coord_lon,
                                                          zen_punct)

                new_file = os.path.join(RM_dir, 'IonRM{hour}.txt'.format(hour=hour))
                utils.write_RM(hour, new_file, B_para, TEC_path, RMS_TEC_path)

                _, _, _, RM_ut, dRM_ut = np.loadtxt(new_file, unpack=True)
                RMs.append(RM_ut)
                dRMs.append(dRM_ut)

            rm_s.append(RMs)
            drm_s.append(dRMs)
            b_para_s.append(B_para)

        self.RMs = np.array(rm_s)
        self.dRMs = np.array(drm_s)
        self.B_paras = np.array(b_para_s)

    def parse_altaz(self):
        '''
        parses ionospheric RM files and assigns values to object

        Returns
        -------
        tuple:
            array[float]: parallel B field array
            array[float]: RM data
            array[float]: RM error data
        '''
        for date in self.times:
            date_str = str(date)
            RM_dir = self.make_rm_dir(date_str)
            RMs = []
            dRMs = []
            for UT in self.UTs:
                new_file = os.path.join(RM_dir, 'IonRM{hour}.txt'.format(hour=hour))
                #B_para differs on a per day basis so any one including last one works
                _, _, B_para, RM_ut, dRM_ut = np.loadtxt(new_file, unpack=True)
                RMs.append(RM_ut)
                dRMs.append(dRM_ut)

            rm_s.append(RMs)
            drm_s.append(dRMs)
            b_para_s.append(B_para)

        self.RMs = np.array(rm_s)
        self.dRMs = np.array(drm_s)
        self.B_paras = np.array(b_para_s)

        return self.B_para, self.RMs, self.dRMs

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

    def map_to_npz(self, time_str, loc_str='PAPER', verbose=False):
        '''
        writes map of particular time to npz file

        Parameters
        ----------
        time_str | str: time and date of map
        loc_str | str: name of obs location to incorporate into filename output
        verbose | Optional[bool]: whether to print values or not
        '''
        #I could fish around in the file read to get npix and then re-loop, but why not just be lazy sometimes
        rng = np.arange(24)
        final_TEC, final_rm, final_drm, ra, dec = np.zeros((rng.shape[0], self.npix)),\
                                                  np.zeros((rng.shape[0], self.npix)),\
                                                  np.zeros((rng.shape[0], self.npix)),\
                                                  np.zeros((rng.shape[0], self.npix)),\
                                                  np.zeros((rng.shape[0], self.npix))
        RM_dir = os.path.join(self.rm_dir, '{date}'.format(date=time_str))
        for UT in rng:
            rm_file = os.path.join(RM_dir, 'IonRM{num}.txt'.format(num=utils.std_hour(UT, verbose=verbose)))
            radec_file = os.path.join(RM_dir, 'radec{num}.txt'.format(num=utils.std_hour(UT, verbose=verbose)))

            _, TEC, B, RM, dRM = np.loadtxt(rm_file, unpack=True)
            RA, DEC = np.loadtxt(radec_file, unpack=True)

            final_TEC[UT, :] = TEC
            final_rm[UT, :] = RM
            final_drm[UT, :] = dRM
            ra[UT,:] = RA
            dec[UT,:] = DEC

        f_name = ''.join((time_str, '_', loc_str, '.npz'))
        npz_dir = os.path.join(rad.root_dir, 'npz')
        if not os.path.exists(npz_dir):
            os.mkdir(npz_dir)
        npz_file = os.path.join(npz_dir, f_name)
        if verbose:
            print('Saving TEC, RM +/- dRM data and RA/Dec mapping to {filename}'.format(filename=npz_file))

        np.savez(npz_file, TEC=final_TEC, RM=final_rm, dRM=final_drm, RA=ra, DEC=dec)

def HERA_RM(date_strs):
    """
    For our convenience: built-in generator for the PAPER/HERA site in the Karoo RQZ, South Africa
    """
    lat_str = '30d43m17.5ss'
    lon_str = '21d25m41.9se'
    height = 1073 # XXX
    return IonoMap(lat_str=lat_str, lon_str=lon_str, date_strs=date_strs, height=height)
