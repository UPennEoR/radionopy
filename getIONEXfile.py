# Generate the necessary IONEX filename for a given day, and fetch it via ftp
from __future__ import print_function
import os
import datetime
import ftplib

# From the little script from ionFR of the same name
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

def getIONEXfile(year, month, day):
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

if __name__ == '__main__':
	#os.system(' '.join('gunzip', IONEX_file_Z))
