import datetime
import numpy as np

# From the little script of the same name in IONEX
# The script is not a function.  
def IONEXFileNeeded(year, month, day):
	time_str = '{year} {month} {day}'.format(year=year, month=month, day=day)
	dayofyear = datetime.datetime.strptime(time_str, '%Y %m %d').timetuple().tm_yday

	if dayofyear < 10:
		dayofyear = '00{dayofyear}'.format(dayofyear=dayofyear)
	elif dayofyear < 100 and dayofyear >= 10:
		dayofyear = '0{dayofyear}'.format(dayofyear=dayofyear)

	# Outputing the name of the IONEX file you require
	ionex_file = 'CODG{dayofyear}0.{yearend}I'.format(dayofyear=dayofyear, yearend=str(year)[2:4])

	return ionex_file

def makeDateString(year, month, day):
	date = '{year}-{month}-{day}T00:00:00'.format(year=year, month=month, day=day)

	return date

def makeIonFRMCommand(RA, Dec, lat, lng, time, ionexfile):
	ionFRMcommand = 'ionFRM.py {RA}{Dec} {lat} {lng} {time} IONEX_Data/ {ionexfile}'.format(RA=RA, Dec=Dec, lat=lat, time=time,
																							lng=lng, ionexfile=ionexfile)

	return ionFRMcommand

def readIonFRtxt():
	UT, TEC, B, RM, dRM = np.loadtxt('IonRM.txt', unpack=True)

	return {'UT': UT,'TEC': TEC,'B': B,'RM': RM,'dRM': dRM}
