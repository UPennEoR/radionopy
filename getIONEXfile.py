# Generate the necessary IONEX filename for a given day, and fetch it via ftp
from __future__ import print_function
import os
import datetime
import ftplib

# From the little script from ionFR of the same name
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

def getIONEXfile(year, month, day):
	server = 'ftp.unibe.ch'

	ftp_dir = os.path.join('aiub/CODE/', year)
	ionexfile = IONEXFileNeeded(year, month, day)
	ionexfileZ = ''.join((ionexfile, '.Z'))

	getting_file_str = 'Retrieving {ionexfileZ} for {day} {month} {year}'.format(ionexfileZ=ionexfileZ, day=day, month=month, year=year)
	print(getting_file_str)

	ftp = ftplib.FTP(server, 'anonymous', 'jaguirre@sas.upenn.edu')
	ftp.cwd(ftp_dir)
	ftp.retrbinary(' '.join(('RETR', ionexfileZ)), open(ionexfileZ, 'wb').write)
	ftp.quit()

	return True

if __name__ == '__main__':
	#os.system(' '.join('gunzip', ionexfileZ))
