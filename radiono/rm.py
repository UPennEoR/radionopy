'''
radiono.rm

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | Module used to create RM object

Classes
-------
RM | allows user to grab data from IONEX files by altaz or radec
'''
from __future__ import print_function
import os
import healpy as hp
import numpy as np
import radiono as rad
from radiono.scripts import ion_altaz as ia, ion_radec as ir

class RM(object):
    '''
    superobject to do altaz or radec RM value generation
    gets the RM based on input

    Methods
    -------
    _radec | gets B_para, RM, and dRM arrays by using ra and dec
    _hp_arr | gets altitude and azimuth arrays for healpix map
    _altaz | gets B_para, RM, and dRM by using alt and az
    '''
    def __init__(self, lat_str, lon_str, time_strs, height=0, ionex_dir=rad.ionex_dir, rm_dir=rad.rm_dir):
        '''
        initalizes RM object

        Parameters
        ----------
        lat_str | str: latitude that the observation was taken at
        lon_str | str: longitude that the observation was taken at
        time_strs | list[str]: list of dates for observation
        height | Optional[float]: height observation taken at in meters
        ionex_dir | Optional[str]: directory in which ionex files are / will be located
        rm_dir | Optional[str]: directory in which RM data files are / will be located
        '''
        self.lat_str = lat_str
        self.lon_str = lon_str
        self.time_strs = time_strs
        self.height = height
        self.ionex_dir = ionex_dir
        self.rm_dir = rm_dir
        self.nside = 16

    def _radec(self, ra_strs, dec_strs, UTs):
        '''
        outputs RM data from radec calculation

        Parameters
        ----------
        ra_strs | list[str]: list of RAs for observation, corresponds by element to dec_strs
        dec_strs | list[str]: list of DECs for observation, corresponds by element for ra_strs
        UTs | array[int]: hours to cycle through

        Returns
        -------
        tuple:
            array[float]: parallel B field array
            array[float]: RM data
            array[float]: RM error data
        '''
        ir.ion_RM(ra_strs, dec_strs,
                  self.lat_str, self.lon_str,
                  self.time_strs, self.height,
                  self.ionex_dir, self.rm_dir)

        rm_s = []
        drm_s = []
        for UT in UTs:
            data_file = os.path.join(self.rm_dir, 'IonRM{hour}.txt'.format(hour=std_hour(UT, verbose=False)))
            _, _, B_para, RM_add, dRM_add = np.loadtxt(data_file, unpack=True)
            rm_s.append(RM_add)
            drm_s.append(dRM_add)

        RMs = np.array(rm_s)
        dRMs = np.array(drm_s)

        return B_para, RMs, dRMs

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

    def _altaz(self, time_str):
        '''
        outputs RM data from altaz calculation

        Parameters
        ----------
        time_str | str: string representation of date and time

        Returns
        -------
        tuple:
            array[float]: parallel B field array
            array[float]: RM data
            array[float]: RM error data
        '''
        alt_src, az_src = self._hp_arr()
        rm_s = []
        drm_s = []
        b_para_s = []
        for time_str in self.time_strs:
            B_para, RM_add, dRM_add = ia.ion_RM(time_str, self.lat_str, self.lon_str, alt_src, az_src, verbose=False)
            rm_s.append(RM_add)
            drm_s.append(dRM_add)
            b_para_s.append(B_para)

        RMs = np.array(rm_s)
        dRMs = np.array(drm_s)
        B_paras = np.array(b_para_s)

        return B_paras, RMs, dRMs

if __name__ == '__main__':
    print('This is not a script anymore')
