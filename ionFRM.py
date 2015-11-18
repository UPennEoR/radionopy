import os
import sys
sys.path.append(os.path.join(base_path, 'ionFR/SiderealPackage'))
sys.path.append(os.path.join(base_path, 'ionFR/PunctureIonosphereCoord'))
sys.path.append(os.path.join(base_path, 'ionFR/IONEX'))
import rdalaz
import ippcoor
import teccalc
import tecrmscalc
import ionheight
import scipy

import subprocess

base_path = '/home/immwa/

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
		zen_s = (scipy.pi / 2.0) - al_s

		# output data only when the altitude of the source is above 0 degrees
		if al_s * (180.0 / scipy.pi) > 0: 

			# Reading the altitude of the Ionosphere in km (from IONEX file)
			alt_ion = ionheight.calc_ion_height(IONEX_name)
			alt_ion = alt_ion * 1000.0 # km to m

			# Alt and AZ coordinates of the Ionospheric piercing point
			# Lon and Lat distances wrt the location of the antenna are also 
			# calculated (radians)
			off_lat, off_lon, az_punct, zen_punct = ippcoor.punc_ion_offset(lat_o, az_s, zen_s, alt_ion)
			al_s_punct = (scipy.pi / 2.0) - zen_punct

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

			TEC_arr = teccalc.calcTEC(lat_val * (lat_o + off_lat) * 180.0 / scipy.pi, lon_val * (lon_o + off_lon) * 180.0 / scipy.pi, IONEX_name)
			RMS_TEC_arr = tecrmscalc.calcRMSTEC(lat_val * (lat_o + off_lat) * 180.0 / scipy.pi, lon_val * (lon_o + off_lon) * 180.0 / pi,
												IONEX_name)

			VTEC = TEC_arr[int(hour)]
			TEC_path = VTEC * TEC2m2 / math.cos(zen_punct) # from vertical TEC to line of sight TEC
			VRMS_TEC = RMS_TEC_arr[int(hour)]
			RMS_TEC_path = VRMS_TEC * TEC2m2 / math.cos(zen_punct) # from vertical RMS TEC to line of sight RMS TEC

			input_file = os.path.join(base_path, 'ionFR/IGRF/geomag70_linux/input.txt')
			output_file = os.path.join(base_path, 'ionFR/IGRF/geomag70_linux/output.txt')

			#uses lat_val, lon_val from above
			# Calculation of the total magnetic field along the line of sight at the IPP
			with open(input_file, 'w') as f:
				f.write('{year},{month},{day} C K{sky_rad} {ipp_lat} {ipp_lon}'.format(year=year, month=month, day=day,
																						sky_rad=(earth_radius + alt_ion) / 1000.0,
																						ipp_lat=lat_val * (lat_o + off_lat) * 180.0 / scipy.pi,
																						ipp_lon=lon_val * (lon_o + off_lon) * 180.0 / scipy.pi)

			#XXX runs the geomag exe script
			script_name = os.path.join(base_path, 'ionFR/IGRF/geomag70_linux/geomag70.exe')
			script_data = os.path.join(base_path, 'ionFR/IGRF/geomag70_linux/IGRF11.COF')
			script_option = 'f'
			subprocess.call([script_name, script_data, script_option, input_file, output_file])

			with open(output_file, 'w') as g:
				data = g.readlines()

			#remove files once used
			os.remove(input_file)
			os.remove(output_file)

			x_field, y_field, z_field = [abs(float(field_data)) * pow(10, -9) * tesla_to_gauss for field_data in data[1].split()[10:13]]
			tot_field = z_field * math.cos(zen_punct) +\
						y_field * math.sin(zen_punct) * math.sin(az_punct) +\
						x_field * math.sin(zen_punct) * math.cos(az_punct)

			# Saving the Ionosheric RM and its corresponding
			# rms value to a file for the given 'hour' value
			IFR = 2.6 * pow(10, -17) * tot_field * TEC_path
			RMS_IFR = 2.6 * pow(10, -17) * tot_field * RMS_TEC_path
			with open(os.path.join(os.getcwd(), 'IonRM.txt'), 'a') as f:
				f.write('{hour} {TEC_path} {tot_field} {IFR} {RMS_IFR}\n'.format(hour=hour, TEC_path=TEC_path, tot_field=tot_field,
																					IFR=IFR, RMS_IFR=RMS_IFR))
