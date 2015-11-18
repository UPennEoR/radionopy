#!/usr/bin/env python
#================================================================
# conjd: Convert Julian date to date and time
#  For documentation, see:
#	 http://www.nmt.edu/tcc/help/lang/python/examples/sidereal/ims/
#----------------------------------------------------------------
#================================================================
# Imports
#----------------------------------------------------------------
from __future__ import print_function
import sys
import sidereal as sdr
#--- main

def main():
	'''
	conjd main program.
	'''

	#-- 1 --
	# [ if sys.argv[1:] is a single float ->
	#	 j := that float
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	arg_list = sys.argv[1:]
	if len(arg_list) != 1:
		usage('Wrong argument count.')
	else:
		try:
			j = float(arg_list[0])
		except ValueError, detail:
			usage('Invalid argument: {detail}'.format(detail=detail))

	#-- 2 --
	# [ jd := a JulianDate instance for Julian date j ]
	jd = sdr.JulianDate(j)

	#-- 3 --
	# [ dt := jd as a datetime.datetime instance ]
	dt = jd.datetime()

	#-- 4 --
	# [ sys.stdout +:= dt in ISO form ]
	print(str(dt))
#--- usage

def usage(*L):
	'''
	Write a usage message and stop.

	 [ L is a list of strings ->
		 sys.stderr +:= (usage message) + (joined elements of L)
		 stop execution ]
	'''
	print(>>sys.stderr, '*** Usage:')
	print(>>sys.stderr, '***  conjd NNNNNNN.NN...')
	print(>>sys.stderr, '*** where NNNNNNN.NN is the Julian date.')
	print(>>sys.stderr, '*** Error: {msg}'.format(msg=''.join(L)))
	raise SystemExit
#================================================================
# Epilogue
#----------------------------------------------------------------

if __name__ == '__main__':
	main()
