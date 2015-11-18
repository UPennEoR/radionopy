#!/usr/bin/env python
#================================================================
# jd: Convert date and time to Julian date
#  For documentation, see:
#	 http://www.nmt.edu/tcc/help/lang/python/examples/sidereal/ims/
#----------------------------------------------------------------
#================================================================
# Imports
#----------------------------------------------------------------
from __future__ import print_function
import sys
import sidereal as sdr
# - - -  m a i n

def main():
	'''
	jd main program.
	'''

	#-- 1 --
	# [ if the arguments in sys.argv are valid ->
	#	 dt := a datetime.datetime instance representing the
	#			 date and time expressed in those arguments ]
	dt = arg_check()

	#-- 2 --
	# [ jd := a JulianDate instance representing dt ]
	jd = sdr.JulianDate.from_datetime(dt)

	#-- 3 --
	print(float(jd))
#--- arg_check

def arg_check():
	'''Check and convert the command line argument(s).
	'''
	#-- 1 --
	# [ arg_list := the command line arguments ]
	arg_list = sys.argv[1:]

	#-- 2 --
	# [ if (len(arg_list)==1) and arg_list[0] is a valid
	#  date-time string ->
	#	 dt := that date-time as a datetime.datetime instance
	#  else if (len(arg_list)==2) and (arg_list[0] is a valid
	#  date) and (arg_list[1] is a valid time) ->
	#	 dt := a datetime.datetime representing that date
	#			 and time
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	if len(arg_list) == 1:
		try:
			dt = sdr.parse_datetime(arg_list[0])
		except SyntaxError, detail:
			usage('Invalid date-time: {detail}'.format(detail=detail))
	elif len(arg_list) == 2:
		try:
			date = sdr.parse_date(arg_list[0])
		except SyntaxError, detail:
			usage('Invalid date: {detail}'.format(detail=detail))
		try:
			time = sdr.parse_time(arg_list[1])
		except SyntaxError, detail:
			usage('Invalid time: {detail}'.format(detail=detail))
		dt = date.combine(date, time)
	else:
		usage('Incorrect number of arguments.')

	#-- 3 --
	return dt
# - - -  u s a g e

def usage(*L):
	'''
	Print a usage message and stop.

	 [ L is a list of strings ->
		 sys.stderr +:= (usage message) + (elements of L,
						  concatenated)
		 stop execution ]
	'''
	print(>>sys.stderr, '*** Usage:')
	print(>>sys.stderr, '***  jd yyyy-mm-dd[Thh[:mm[:ss]]]')
	print(>>sys.stderr, '*** or:')
	print(>>sys.stderr, '***  jd yyyy-mm-dd hh[:mm[:ss]]')
	print(>>sys.stderr, '*** Error: {msg}'.format(msg=''.join(L)))
	raise SystemExit
#================================================================
# Epilogue
#----------------------------------------------------------------

if __name__ == '__main__':
	main()
