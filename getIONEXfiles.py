# There's always a beter way, but this will get the files I need to
# reply to the referee

# 2455903 until 2455985
# 7 Dec 2011 to 27 Feb 2012
import ftplib
import datetime
import os

# From the little script of the same name
def IONEXFileNeeded(year,month,day):
    dayofyear = datetime.datetime.strptime(''+str(year)+' '+str(month)+' '+str(day)+'', '%Y %m %d').timetuple().tm_yday

    if dayofyear < 10:
        dayofyear = '00'+str(dayofyear)
    if dayofyear < 100 and dayofyear >= 10:
        dayofyear = '0'+str(dayofyear)

    # Outputing the name of the IONEX file you require
    file =  'CODG'+str(dayofyear)+'0.'+str(list(str(year))[2])+''+str(list(str(year))[3])+'I'
    return file

year = '2012'
month = '01'
# pad it with one day on each side
dec = range(6,32)
jan = range(1,32)
feb = range(1,29)

day = []
month = []
year = []

for i in dec:
    year.append('2011')
    month.append('12')
    day.append('%.2d'%(i))
for i in jan:
    year.append('2012')
    month.append('01')
    day.append('%.2d'%(i))
for i in feb:
    year.append('2012')
    month.append('02')
    day.append('%.2d'%(i))

server = 'ftp.unibe.ch'

for i in range(len(year)):
    ftp_dir = 'aiub/CODE/'+year[i]

    ionexfile = IONEXFileNeeded(year[i],month[i],day[i])
    ionexfileZ = ionexfile+'.Z'

    print 'Retrieving '+ionexfileZ+' for '+day[i]+' '+month[i]+' '+year[i]
    ftp = ftplib.FTP(server, 'anonymous', 'jaguirre@sas.upenn.edu')
    ftp.cwd(ftp_dir)
    ftp.retrbinary('RETR '+ionexfileZ, open(ionexfileZ, 'wb').write)
    ftp.quit()

#    os.system('gunzip '+ionexfileZ)

