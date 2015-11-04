# Generate the necessary IONEX filename for a given day, and fetch it via ftp

import ftplib
import datetime
import os

# From the little script from ionFR of the same name
def IONEXFileNeeded(year,month,day):
    dayofyear = datetime.datetime.strptime(''+str(year)+' '+str(month)+' '+str(day)+'', '%Y %m %d').timetuple().tm_yday

    if dayofyear < 10:
        dayofyear = '00'+str(dayofyear)
    if dayofyear < 100 and dayofyear >= 10:
        dayofyear = '0'+str(dayofyear)

    # Outputing the name of the IONEX file you require
    file =  'CODG'+str(dayofyear)+'0.'+str(list(str(year))[2])+''+str(list(str(year))[3])+'I'
    return file

server = 'ftp.unibe.ch'

def getIONEXfile(year,month,day):

    ftp_dir = 'aiub/CODE/'+year
    ionexfile = IONEXFileNeeded(year,month,day)
    ionexfileZ = ionexfile+'.Z'

    print 'Retrieving '+ionexfileZ+' for '+day+' '+month+' '+year
    ftp = ftplib.FTP(server, 'anonymous', 'jaguirre@sas.upenn.edu')
    ftp.cwd(ftp_dir)
    ftp.retrbinary('RETR '+ionexfileZ, open(ionexfileZ, 'wb').write)
    ftp.quit()

    return True

#    os.system('gunzip '+ionexfileZ)

