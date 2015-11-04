import datetime
import numpy as np

# From the little script of the same name in IONEX
# The script is not a function.  
def IONEXFileNeeded(year,month,day):
    dayofyear = datetime.datetime.strptime(''+str(year)+' '+str(month)+' '+str(day)+'', '%Y %m %d').timetuple().tm_yday

    if dayofyear < 10:
        dayofyear = '00'+str(dayofyear)
    if dayofyear < 100 and dayofyear >= 10:
        dayofyear = '0'+str(dayofyear)

    # Outputing the name of the IONEX file you require
    file =  'CODG'+str(dayofyear)+'0.'+str(list(str(year))[2])+''+str(list(str(year))[3])+'I'
    return file

def makeDateString(year,month,day):
    date = str(year)+'-'+str(month)+'-'+str(day)+'T00:00:00'
    return date

def makeIonFRMCommand(RA,Dec,lat,lng,time,ionexfile):
    ionFRMcommand = 'ionFRM.py '+RA+Dec+' '+lat+' '+lng+' '+time+' '+'IONEX_Data/'+ionexfile
    return ionFRMcommand

def readIonFRtxt():
    UT,TEC,B,RM,dRM = np.loadtxt('IonRM.txt',unpack=True)
    return {'UT':UT,'TEC':TEC,'B':B,'RM':RM,'dRM':dRM}
