import os
import sys
sys.path.append(os.path.join(base_path, 'ionFR/SiderealPackage'))
sys.path.append(os.path.join(base_path, 'ionFR/PunctureIonosphereCoord'))
sys.path.append(os.path.join(base_path, 'ionFR/IONEX'))
import rdalaz
import teccalc
import tecrmscalc
import numpy as np
import subprocess

base_path = os.path.expanduser('~')

def IONEX_file_needed(year, month, day):
	time_str = '{year} {month} {day}'.format(year=year, month=month, day=day)
	day_of_year = datetime.datetime.strptime(time_str, '%Y %m %d').timetuple().tm_yday

	if day_of_year < 10:
		day_of_year = '00{day_of_year}'.format(day_of_year=day_of_year)
	elif 10 <= day_of_year < 100:
		day_of_year = '0{day_of_year}'.format(day_of_year=day_of_year)

	# Outputing the name of the IONEX file you require
	ionex_file = 'CODG{day_of_year}0.{year_end}I'.format(day_of_year=day_of_year, year_end=str(year)[2:4])

	return ionex_file

def get_IONEX_file(year, month, day):
	server = 'ftp.unibe.ch'

	ftp_dir = os.path.join('aiub/CODE/', year)
	IONEX_file = IONEX_file_needed(year, month, day)
	IONEX_file_Z = ''.join((IONEX_file, '.Z'))

	getting_file_str = 'Retrieving {IONEX_file_Z} for {day} {month} {year}'.format(IONEX_file_Z=IONEX_file_Z, day=day, month=month, year=year)
	print(getting_file_str)

	ftp = ftplib.FTP(server, 'anonymous', 'jaguirre@sas.upenn.edu')
	ftp.cwd(ftp_dir)
	ftp.retrbinary(' '.join(('RETR', IONEX_file_Z)), open(IONEX_file_Z, 'wb').write)
	ftp.quit()

	return True

def read_IONEX_TEC(filename):
	#==========================================================================
	# Reading and storing only the TEC values of 1 day
	# (13 maps) into a 3D array

	# Opening and reading the IONEX file into memory
	with open(filename, 'r') as read_file:
		linestring = read_file.read()
		IONEX_list = linestring.split('\n')

	# creating a new array without the header and only
	# with the TEC maps
	add = 0 
	new_IONEX_list = []
	for file_data in IONEX_list:
		if file_data.split()[-2:] == ['RMS', 'MAP']:
			add = 0
		elif file_data.split()[-2:] == ['IN', 'FILE']:
			number_of_maps = float(file_data.split()[0])

		if add == 1:
			new_IONEX_list.append(file_data)

		if file_data.split()[0] == 'END' and file_data.split()[2] == 'HEADER':
			add = 1

		if file_data.split()[-1] == 'DHGT':
			ion_h = float(file_data.split()[0])
		elif file_data.split()[-1] == 'DLAT':
			start_lat, end_lat, step_lat = [float(data_item) for data_item in file_data[:2]]
		elif file_data.split()[-1] == 'DLON':
			start_lon, end_lon, step_lon = [float(data_item) for data_item in file_data[:2]]

	# Variables that indicate the number of points in Lat. and Lon.
	points_lon = ((end_lon - start_lon) / step_lon) + 1
	points_lat = ((end_lat - start_lat) / step_lat) + 1

	print(start_lon, end_lon, step_lon)
	print(start_lat, end_lat, step_lat)
	print(points_lon, points_lat)

	# What are the Lat/Lon coords?
	longitude = np.linspace(start_lon, end_lon, num=points_lon)
	latitude = np.linspace(start_lat, end_lat, num=points_lat)

	# 3D array that will contain TEC values only
	a = np.zeros((number_of_maps, points_lat, points_lon))

	# Selecting only the TEC values to store in the 3-D array
	counter_maps = 1
	for i in range(len(new_IONEX_list)):
		# Pointing to first map (out of 13 maps)
		# then by changing 'counter_maps' the other
		# maps are selected
		if new_IONEX_list[i].split()[0] == str(counter_maps) and new_IONEX_list[i].split()[-4] == 'START':
			# pointing the starting latitude
			# then by changing 'counter_lat' we select
			# TEC data at other latitudes within
			# the selected map
			counter_lat = 0
			new_start_lat = float(str(start_lat))
			for item_lat in range(int(points_lat)):
				if new_IONEX_list[i + 2 + counter_lat].split()[0].split('-')[0] == str(new_start_lat)\
				or '-' + new_IONEX_list[i + 2 + counter_lat].split()[0].split('-')[1] == str(new_start_lat):
					# Adding to array 'a' a line of latitude TEC data
					# we account for the TEC values at negative latitudes
					counter_lon = 0
					for count_num in range(3, 8):
						list_index = i + count_num + counter_lat
						for new_IONEX_item in new_IONEX_list[list_index].split():
							a[counter_maps - 1, item_lat, counter_lon] = new_IONEX_item
							counter_lon = counter_lon + 1
				counter_lat = counter_lat + 6
				new_start_lat = new_start_lat + step_lat
			counter_maps = counter_maps + 1

	return {'TEC': np.array(a), 'lat': latitude, 'lon': longitude}
    return {'TEC': np.array(a), 'lat': latitude, 'lon': longitude, 'AltIon': ion_h * 1000.0}
	#==========================================================================

def punc_ion_offset(lat_obs, az_sou, ze_sou, alt_ion):
	radius_earth = 6371000.0 # in meters

	# The 2-D sine rule gives the zenith angle at the
	# Ionospheric piercing point
	zen_punc = np.asin((radius_earth * np.sin(ze_sou)) / (radius_earth + alt_ion)) 

	# Use the sum of the internal angles of a triange to determine theta
	theta = ze_sou - zen_punc

	# The cosine rule for spherical triangles gives us the latitude
	# at the IPP
	lat_ion = np.asin(np.sin(lat_obs) * np.cos(theta) + np.cos(lat_obs) * np.sin(theta) * np.cos(az_sou)) 
	d_lat = lat_ion - lat_obs # latitude difference

	# Longitude difference using the 3-D sine rule (or for spherical triangles)
	d_lon = np.asin(np.sin(az_sou) * np.sin(theta) / np.cos(lat_ion))

	# Azimuth at the IPP using the 3-D sine rule
	s_az_ion = np.sin(az_sou) * np.cos(lat_obs) / np.cos(lat_ion)
	az_punc = np.asin(s_az_ion)

	return d_lat, d_lon, az_punc, zen_punc

def calc_ion_height(filename): 
	# opening and reading the IONEX file into memory
	with open(filename, 'r') as read_file:
		linestring = read_file.read()
		IONEX_list = linestring.split('\n')

	for file_data in IONEX_list:
		if file_data.split()[-1] == 'DHGT'
			ion_h = float(file_data.split()[0])

	return ion_h

def B_IGRF(year, month, day, alt_ion, lon_o, lat_o, off_lat, az_punct, zen_punct):
	al_s_punct = (np.pi / 2.0) - zen_punct

	# Calculation of TEC path value for the indicated 'hour' and therefore 
	# at the IPP
	# Calculation of RMS TEC path value (same as the step above)
	if raw_latitude[-1] == 's':
		lat_val = -1
	elif raw_latitude[-1] == 'n':
		lat_val = 1
	if raw_longitude[-1] == 'e':
		lon_val = 1
	elif raw_longitude[-1] == 'w':
		lon_val = -1

	TEC_arr = teccalc.calcTEC(lat_val * (lat_o + off_lat) * 180.0 / np.pi, lon_val * (lon_o + off_lon) * 180.0 / np.pi, IONEX_name)
	RMS_TEC_arr = tecrmscalc.calcRMSTEC(lat_val * (lat_o + off_lat) * 180.0 / np.pi, lon_val * (lon_o + off_lon) * 180.0 / pi,
										IONEX_name)

	VTEC = TEC_arr[int(hour)]
	TEC_path = VTEC * TEC2m2 / np.cos(zen_punct) # from vertical TEC to line of sight TEC
	VRMS_TEC = RMS_TEC_arr[int(hour)]
	RMS_TEC_path = VRMS_TEC * TEC2m2 / np.cos(zen_punct) # from vertical RMS TEC to line of sight RMS TEC

	input_file = os.path.join(base_path, 'ionFR/IGRF/geomag70_linux/input.txt')
	output_file = os.path.join(base_path, 'ionFR/IGRF/geomag70_linux/output.txt')

	#uses lat_val, lon_val from above
	# Calculation of the total magnetic field along the line of sight at the IPP
	with open(input_file, 'w') as f:
		f.write('{year},{month},{day} C K{sky_rad} {ipp_lat} {ipp_lon}'.format(year=year, month=month, day=day,
																				sky_rad=(earth_radius + alt_ion) / 1000.0,
																				ipp_lat=lat_val * (lat_o + off_lat) * 180.0 / np.pi,
																				ipp_lon=lon_val * (lon_o + off_lon) * 180.0 / np.pi)

	#XXX runs the geomag exe script
	script_name = os.path.join(base_path, 'ionFR/IGRF/geomag70_linux/geomag70.exe')
	script_data = os.path.join(base_path, 'ionFR/IGRF/geomag70_linux/IGRF11.COF')
	script_option = 'f'
	subprocess.call([script_name, script_data, script_option, input_file, output_file])

	with open(output_file, 'w') as g:
		data = g.readlines()


		x_field, y_field, z_field = [abs(float(field_data)) * pow(10, -9) * tesla_to_gauss for field_data in data[1].split()[10:13]]
		tot_field = z_field * np.cos(zen_punct) +\
					y_field * np.sin(zen_punct) * np.sin(az_punct) +\
					x_field * np.sin(zen_punct) * np.cos(az_punct)

		#remove files once used
		os.remove(input_file)
		os.remove(output_file)

	return tot_field

if __name__ == '__main__':
	# Defining some variables for further use
	TECU = pow(10, 16)
	TEC2m2 = 0.1 * TECU
	earth_radius = 6371000.0 # in meters
	tesla_to_gauss = pow(10, 4)

	# Cheking the arguments are given correctly
	arg_list  =  sys.argv[1:]
	if len(arg_list) != 5:
		usage('Incorrect command line argument count.')
	else:
		raw_ra_dec, raw_latitude, raw_longitude, raw_d_time, IONEX_name = arg_list

	## Nominally try to reproduce the output of this command
	## ionFRM.py 16h50m04.0s+79d11m25.0s 52d54m54.64sn 6d36m16.04se 2004-05-19T00:00:00 CODG1400.04I
	## Echo back what he has ... 
	#RAstr = '16h50m04.0s'
	#DECstr = '+79d11m25.0s'
	#LatStr = '52d54m54.64sn'
	#LonStr = '6d36m16.04se'
	#TimeStr = '2004-05-19T00:00:00' # This will actually work as input to the astropy Time function
	#IONEXfile = 'CODG1400.04I'

	# predict the ionospheric RM for every hour within a day 
	for h in range(24):
		if h < 10:
			raw_time = '{date_time}T0{hour}:00:00'.format(date_time=raw_d_time.split('T')[0], hour=h)
		else:
			raw_time = '{date_time}T{hour}:00:00'.format(date_time=raw_d_time.split('T')[0], hour=h)
		
		hour = raw_time.split('T')[1].split(':')[0]
		date = raw_time.split('T')[0].split('-')
		year, month, day = raw_time.split('T')[0].split('-')[:2]

		# RA and Dec (of the source in degrees) to Alt and Az (radians)
		az_s, al_s, ha, lat_o, lon_o = rdalaz.alt_az(raw_time)
		zen_s = (np.pi / 2.0) - al_s

		# output data only when the altitude of the source is above 0 degrees
		if al_s * (180.0 / np.pi) > 0: 

			# Reading the altitude of the Ionosphere in km (from IONEX file)
			alt_ion = calc_ion_height(IONEX_name)
			alt_ion = alt_ion * 1000.0 # km to m

			# Alt and AZ coordinates of the Ionospheric piercing point
			# Lon and Lat distances wrt the location of the antenna are also 
			# calculated (radians)
			off_lat, off_lon, az_punct, zen_punct = punc_ion_offset(lat_o, az_s, zen_s, alt_ion)
			tot_field = B_IGRF(year, month, day, alt_ion, lon_o, lat_o, off_lat, az_punct, zen_punct)

			# Saving the Ionosheric RM and its corresponding
			# rms value to a file for the given 'hour' value
			IFR = 2.6 * pow(10, -17) * tot_field * TEC_path
			RMS_IFR = 2.6 * pow(10, -17) * tot_field * RMS_TEC_path
			with open(os.path.join(os.getcwd(), 'IonRM.txt'), 'a') as f:
				f.write('{hour} {TEC_path} {tot_field} {IFR} {RMS_IFR}\n'.format(hour=hour, TEC_path=TEC_path, tot_field=tot_field,
																					IFR=IFR, RMS_IFR=RMS_IFR))
