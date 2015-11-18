import datetime
import numpy as np

# From the little script of the same name in IONEX
# The script is not a function.  
def IONEX_file_needed(year, month, day):
	time_str = '{year} {month} {day}'.format(year=year, month=month, day=day)
	day_of_year = datetime.datetime.strptime(time_str, '%Y %m %d').timetuple().tm_yday

	if day_of_year < 10:
		day_of_year = '00{day_of_year}'.format(day_of_year=day_of_year)
	elif day_of_year < 100 and day_of_year >= 10:
		day_of_year = '0{day_of_year}'.format(day_of_year=day_of_year)

	# Outputing the name of the IONEX file you require
	ionex_file = 'CODG{day_of_year}0.{year_end}I'.format(day_of_year=day_of_year, year_end=str(year)[2:4])

	return ionex_file

def make_date_string(year, month, day):
	date = '{year}-{month}-{day}T00:00:00'.format(year=year, month=month, day=day)

	return date

def make_ion_frm_command(ra, dec, lat, lng, time, ionexfile):
	ionFRMcommand = 'ionFRM.py {ra}{dec} {lat} {lng} {time} IONEX_Data/ {ionexfile}'.format(ra=ra, dec=dec, lat=lat, time=time,
																							lng=lng, ionexfile=ionexfile)

	return ionFRMcommand

def read_ion_fr():
	UT, TEC, B, RM, dRM = np.loadtxt('IonRM.txt', unpack=True)

	return {'UT': UT,'TEC': TEC,'B': B,'RM': RM,'dRM': dRM}
