'''sidereal.py: A Python module for astronomical calculations.

 For documentation, see:
	http://www.nmt.edu/tcc/help/lang/python/examples/sidereal/ims/
'''
#================================================================
# Imports
#----------------------------------------------------------------
from __future__ import print_function
import math
import re
import datetime
#================================================================
# Manifest constants
#----------------------------------------------------------------

FIRST_GREGORIAN_YEAR = 1583
TWO_PI = 2.0 * math.pi
PI_OVER_12 = math.pi / 12.0
JULIAN_BIAS = 2200000	# 2,200,000
SIDEREAL_A = 0.0657098
FLOAT_PAT = re.compile(r'\d+'		 # Matches one or more digits
						r'('			# Start optional fraction
						r'[.]'		 # Matches the decimal point
						r'\d+'		 # Matches one or more digits
						r')?')		 # End optional group
D_PAT = re.compile(r'[dD]')
M_PAT = re.compile(r'[mM]')
S_PAT = re.compile(r'[sS]')
H_PAT = re.compile(r'[hH]')
NS_PAT = re.compile(r'[nNsS]')
EW_PAT = re.compile(r'[eEwW]')
# - - -  h o u r s T o R a d i a n s

def hours_to_radians(hours):
	'''Convert hours (15 degrees) to radians.
	'''
	return hours * PI_OVER_12
# - - -  r a d i a n s T o H o u r s

def radians_to_hours(radians):
	'''Convert radians to hours (15 degrees).
	'''
	return radians / PI_OVER_12
# - - -  h o u r A n g l e T o R A

def hour_angle_to_ra(h, ut, e_lon):
	'''Convert hour angle to right ascension.

	 [ (h is an hour angle in radians as a float) and
		(ut is a timestamp as a datetime.datetime instance) and
		(e_lon is an east longitude in radians) ->
		 return the right ascension in radians corresponding
		 to that hour angle at that time and location ]
	'''
	#-- 1 --
	# [ gst := the Greenwich Sidereal Time equivalent to
	#			ut, as a SiderealTime instance ]
	gst = SiderealTime.from_datetime(ut)
	#-- 2 --
	# [ lst := the local time corresponding to gst at
	#			longitude e_lon ]
	lst = gst.lst(e_lon)
	#-- 3 --
	# [ alpha := lst - h, normalized to [0, 2 * math.pi) ]
	alpha = (lst.radians - h) % TWO_PI

	#-- 4 --
	return alpha
# - - -  r a T o H o u r A n g l e

def ra_to_hour_angle(ra, ut, e_lon):
	'''Convert right ascension to hour angle.

	 [ (ra is a right ascension in radians as a float) and
		(ut is a timestamp as a datetime.datetime instance) and
		(e_lon is an east longitude in radians) ->
		 return the hour angle in radians at that time and
		 location corresponding to that right ascension ]
	'''
	#-- 1 --
	# [ gst := the Greenwich Sidereal Time equivalent to
	#			ut, as a SiderealTime instance ]
	gst = SiderealTime.from_datetime(ut)

	#-- 2 --
	# [ lst := the local time corresponding to gst at
	#			longitude e_lon ]
	lst = gst.lst(e_lon)
	#-- 3 --
	# [ h := lst - ra, normalized to [0, 2 * math.pi) ]
	h = (lst.radians - ra) % TWO_PI

	#-- 4 --
	return h
# - - -  d a y N o

def day_num(dt):
	'''Compute the day number within the year.

	 [ dt is a date as a datetime.datetime or datetime.date ->
		 return the number of days between dt and Dec. 31 of
		 the preceding year ]
	'''
	#-- 1 --
	# [ date_ord := proleptic Gregorian ordinal of dt
	#  jan_1_ord := proleptic Gregorian ordinal of January 1
	#				of year (dt.year) ]
	date_ord = dt.toordinal()
	jan_1_ord = datetime.date(dt.year, 1, 1).toordinal()

	#-- 2 --
	return date_ord - jan_1_ord + 1
# - - -  p a r s e D a t e t i m e

T_PATTERN = re.compile('[tT]')

def parse_datetime(s):
	'''Parse a date with optional time.

	 [ s is a string ->
		 if s is a valid date with optional time ->
			return that timestamp as a datetime.datetime instance
		 else -> raise SyntaxError ]
	'''
	#-- 1 --
	# [ if s contains 'T' or 't' ->
	#	 raw_date := s up to the first such character
	#	 raw_time := s from just after the first such
	#				 character to the end
	#  else ->
	#	 raw_date := s
	#	 raw_time := None ]
	m = T_PATTERN.search(s)
	if m is None:
		raw_date = s
		raw_time = None
	else:
		raw_date = s[:m.start()]
		raw_time = s[m.end():]
	#-- 2 --
	# [ if raw_date is a valid date ->
	#	 date_part := raw_date as a datetime.datetime instance
	#  else -> raise SyntaxError ]
	date_part = parse_date(raw_date)
	#-- 3 --
	# [ if raw_time is None ->
	#	 time_part := 00:00 as a datetime.time
	#  else if raw_time is valid ->
	#	 time_part := raw_time as a datetime.time
	#  else -> raise SyntaxError ]
	if raw_time is None:
		time_part = datetime.time(0, 0)
	else:
		time_part = parse_time(raw_time)
	#-- 4 --
	return datetime.datetime.combine(date_part, time_part)
# - - -  p a r s e D a t e

YEAR_FIELD = 'Y'
MONTH_FIELD = 'M'
DAY_FIELD = 'D'

date_re = (r'('			# Begin YEAR_FIELD
			r'?P<%s>'	  # Name this group YEAR_FIELD
			r'\d{4}'		# Match exactly four digits
			r')'			# End YEAR_FIELD
			r'\-'		  # Matches one hyphen
			r'('			# Begin MONTH_FIELD
			r'?P<%s>'	  # Name this group MONTH_FIELD
			r'\d{1,2}'	 # Matches one or two digits
			r')'			# End MONTH_FIELD
			r'\-'		  # Matches '-'
			r'('			# Begin DAY_FIELD
			r'?P<%s>'	  # Name this group DAY_FIELD
			r'\d{1,2}'	 # Matches one or two digits
			r')'			# End DAY_FIELD
			r'$'			# Make sure all characters match
			) % (YEAR_FIELD, MONTH_FIELD, DAY_FIELD)
DATE_PAT = re.compile(date_re)

def parse_date(s):
	'''Validate and convert a date in external form.

	 [ s is a string ->
		 if s is a valid external date string ->
			return that date as a datetime.date instance
		 else -> raise SyntaxError ]
	'''
	#-- 1 --
	# [ if DATE_PAT matches s ->
	#	 m := a match instance describing the match
	#  else -> raise SyntaxError ]
	m = DATE_PAT.match(s)
	if m is None:
		raise SyntaxError, ("Date does not have pattern YYYY-DD-MM: '{date_str}'".format(date_str=date_str)
	#-- 2 --
	year  = int(m.group(YEAR_FIELD))
	month = int(m.group(MONTH_FIELD))
	day	= int(m.group(DAY_FIELD))

	#-- 3 --
	return datetime.date(year, month, day)
# - - -  p a r s e T i m e

def parse_time(s):
	'''Validate and convert a time and optional zone.

	 [ s is a string ->
		 if s is a valid time with optional zone suffix ->
			return that time as a datetime.time
		 else -> raise SyntaxError ]
	'''
	#-- 1 -
	# [ if s starts with FLOAT_PAT ->
	#	 dec_hour := matching part of s as a float
	#	 min_tail := part s past the match
	#  else -> raise SyntaxError ]
	dec_hour, min_tail = parse_float(s, 'Hour number')
	#-- 2 --
	# [ if min_tail starts with ':' followed by FLOAT_PAT ->
	#	 dec_min := part matching FLOAT_PAT as a float
	#	 sec_tail := part of min_tail after the match
	#  else if min_tail starts with ':' not followed by
	#  FLOAT_PAT ->
	#	 raise SyntaxError
	#  else ->
	#	 dec_min := 0.0
	#	 sec_tail := min_tail ]
	if min_tail.startswith(':'):
		m = FLOAT_PAT.match(min_tail[1:])
		if m is None:
			raise SyntaxError, ("Expecting minutes: '{min_tail'".format(min_tail=min_tail))
		else:
			dec_min = float(m.group())
			sec_tail = min_tail[m.end()+1:]
	else:
		dec_min = 0.0
		sec_tail = min_tail
	#-- 3 --
	# [ if sec_tail starts with ':' followed by FLOAT_PAT ->
	#	 dec_sec := part matching FLOAT_PAT as a float
	#	 zone_tail := part of sec_tail after the match
	#  else if sec_tail starts with ':' not followed by
	#  FLOAT_PAT ->
	#	 raise SyntaxError
	#  else ->
	#	 dec_sec := 0.0
	#	 zone_tail := sec_tail ]
	if sec_tail.startswith(':'):
		m = FLOAT_PAT.match(sec_tail[1:])
		if m is None:
			raise SyntaxError, ("Expecting seconds: '{sec_tail}'".format(sec_tail=sec_tail))
		else:
			dec_sec = float(m.group())
			zone_tail = sec_tail[m.end()+1:]
	else:
		dec_sec = 0.0
		zone_tail = sec_tail
	#-- 4 --
	# [ if zone_tail is empty ->
	#	 tz := None
	#  else if zone_tail is a valid zone suffix ->
	#	 tz := that zone information as an instance of a class
	#			 that inherits from datetime.tzinfo
	#  else -> raise SyntaxError ]
	if len(zone_tail) == 0:
		tz = None
	else:
		tz = parse_zone(zone_tail)
	#-- 5 --
	# [ hours := dec_hour + dec_min/60.0 + dec_sec/3600.0 ]
	hours = dms_units.mix_to_single((dec_hour, dec_min, dec_sec))
	#-- 6 --
	# [ return a datetime.time representing hours ]
	hh, mm, secs = dms_units.single_to_mix(hours)
	whole_secs, frac_secs = divmod(secs, 1.0)
	ss = int(whole_secs)
	usec = int(frac_secs * 1e6)
	return datetime.time(hh, mm, ss, usec, tz)
# - - -  p a r s e Z o n e

def parse_zone(s):
	'''Validate and convert a time zone suffix.

	 [ s is a string ->
		 if s is a valid time zone suffix ->
			return that zone's information as an instance of
			a class that inherits from datetime.tzinfo
		 else -> raise SyntaxError ]
	'''
	#-- 1 --
	# [ if s starts with '+' or '-' and is a valid fixed-offset
	#  time zone suffix ->
	#	 return that zone's information as a datetime.tzinfo instance
	#  else if is starts with '+' or '-' but is not a valid
	#  fixed-offset time zone suffix ->
	#	 raise SyntaxError
	#  else -> I ]
	if s.startswith('+') or s.startswith('-'):
		return parse_fixed_zone(s)

	#-- 2 --
	# [ if s.upper() is a key in zone_code_map ->
	#	 return the corresponding value
	#  else -> raise SyntaxError ]
	try:
		tz = zone_code_map[s.upper()]
		return tz
	except KeyError:
		raise SyntaxError, ("Unknown time zone code: '{zone_code}'".format(zone_code=s))
# - - -  p a r s e F i x e d Z o n e

HHMM_PAT = re.compile(r'\d{4}'	# Matches exactly four digits
						r'$')		# Be sure everything is matched

def parse_fixed_zone(s):
	'''Convert a +hhmm or -hhmm zone suffix.

	 [ s is a string ->
		 if s is a time zone suffix of the form '+hhmm' or '-hhmm' ->
			return that zone information as an instance of a class
			that inherits from datetime.tzinfo
		 else -> raise SyntaxError ]
	'''
	#-- 1 --
	if s.startswith('+'):
		sign = 1
	elif s.startswith('-'):
		sign = -1
	else:
		raise SyntaxError, ("Expecting zone modifier as {zone_suffix}hhmm: '{time_zone}'".format(zone_suffix=s[0], time_zone=s))

	#-- 2 --
	# [ if s[1:] matches HHMM_PAT ->
	#	 hours := the HH part as an int
	#	 minutes := the MM part as an int
	#  else -> raise SyntaxError ]
	raw_hhmm = s[1:]
	m = HHMM_PAT.match(raw_hhmm)
	if m is None:
		raise SyntaxError, ("Expecting zone modifier as {zone_suffix}hhmm: '{time_zone}'".format(zone_suffix=s[0], time_zone=s))
	else:
		hours = int(raw_hhmm[:2])
		minutes = int(raw_hhmm[2:])

	#-- 3 --
	return FixedZone(sign * hours, sign * minutes, s)

# - - - - -  c l a s s  F i x e d Z o n e

DELTA_ZERO = datetime.timedelta(0)
DELTA_HOUR = datetime.timedelta(hours=1)

class FixedZone(datetime.tzinfo):
	'''Represents a time zone with a fixed offset east of UTC.

	 Exports:
		FixedZone(hours, minutes, name):
		 [ (hours is a signed offset in hours as an int) and
			(minutes is a signed offset in minutes as an int) ->
			 return a new FixedZone instance representing
			 those offsets east of UTC ]
	 State/Invariants:
		.__offset:
		 [ a datetime.timedelta representing self's offset
			east of UTC ]
		.__name:
		 [ as passed to the constructor's name argument ]
	'''
	def __init__(self, hh, mm, name):
		'''Constructor for FixedZone.
		'''
		self.__offset = datetime.timedelta(hours=hh, minutes=mm)
		self.__name = name
	def utc_offset(self, dt):
		'''Return self's offset east of UTC.
		'''
		return self.__offset
	def tzname(self, dt):
		'''Return self's name.
		'''
		return self.__name
	def dst(self, dt):
		'''Return self's daylight time offset.
		'''
		return DELTA_ZERO

def first_sunday_on_or_after(dt):
	'''Find the first Sunday on or after a given date.

	 [ dt is a datetime.date ->
		 return a datetime.date representing the first Sunday
		 on or after dt ]
	'''
	days_to_go = dt.weekday()
	if days_to_go:
		dt += datetime.timedelta(days_to_go)
	return dt

# - - - - -  c l a s s  U S T i m e Z o n e

class USTimeZone(datetime.tzinfo):
	'''Represents a U.S. time zone, with automatic daylight time.

	 Exports:
		USTimeZone(hh, mm, name, std_name, dst_name):
		 [ (hh is an offset east of UTC in hours) and
			(mm is an offset east of UTC in minutes) and
			(name is the composite zone name) and
			(std_name is the non-DST name) and
			(dst_name is the DST name) ->
			 return a new USTimeZone instance with those values ]

	 State/Invariants:
		.__offset:
		 [ self's offset east of UTC as a datetime.timedelta ]
		.__name:	 [ as passed to constructor's name ]
		.__std_name:  [ as passed to constructor's std_name ]
		.__dst_name:  [ as passed to constructor's dst_name ]
	'''
	DST_START_OLD = datetime.datetime(1, 4, 1, 2)
	DST_END_OLD = datetime.datetime(1, 10, 25, 2)
	DST_START_2007 = datetime.datetime(1, 3, 8, 2)
	DST_END_2007 = datetime.datetime(1, 11, 1, 2)
	def __init__(self, hh, mm, name, std_name, dst_name):
		self.__offset = datetime.timedelta(hours=hh, minutes=mm)
		self.__name	 = name
		self.__std_name = std_name
		self.__dst_name = dst_name
	def tzname(self, dt):
		if self.dst(dt):
			return self.__dst_name
		else:
			return self.__std_name
	def utc_offset(self, dt):
		return self.__offset + self.dst(dt)
	def dst(self, dt):
		'''Return the current DST offset.

		 [ dt is a datetime.date ->
			 if daylight time is in effect in self's zone on
			 date dt ->
				return +1 hour as a datetime.timedelta
			 else ->
				return 0 as a datetime.delta ]
		'''
		#-- 1 --
		# [ dt_start := Sunday when DST starts in year dt.year
		#  dt_end	:= Sunday when DST ends in year dt.year ]
		if dt.year >= 2007:
			start_date = self.DST_START_2007.replace(year=dt.year)
			end_date = self.DST_END_2007.replace(year=dt.year)
		else:
			start_date = self.DST_START_OLD.replace(year=dt.year)
			end_date = self.DST_END_OLD.replace(year=dt.year)
		dt_start = first_sunday_on_or_after(start_date)
		dt_end = first_sunday_on_or_after(end_date)
		#-- 2 --
		# [ naive_date := dt with its tzinfo member set to None ]
		naive_date = dt.replace(tzinfo=None)
		#-- 3 --
		# [ if naive_date is in the interval (dt_start, dt_end) ->
		#	 return DELTA_HOUR
		#  else ->
		#	 return DELTA_ZERO ]
		if dt_start <= naive_date < dt_end:
			return DELTA_HOUR
		else:
			return DELTA_ZERO

utc_zone = FixedZone(0, 0, 'UTC')

est_zone = FixedZone(-5, 0, 'EST')
edt_zone = FixedZone(-4, 0, 'EDT')
et_zone  = USTimeZone(-5, 0, 'ET', 'EST', 'EDT')

cst_zone = FixedZone(-6, 0, 'CST')
cdt_zone = FixedZone(-5, 0, 'CDT')
ct_zone  = USTimeZone(-6, 0, 'CT', 'CST', 'CDT')

mst_zone = FixedZone(-7, 0, 'MST')
mdt_zone = FixedZone(-6, 0, 'MDT')
mt_zone  = USTimeZone(-7, 0, 'MT', 'MST', 'MDT')

pst_zone = FixedZone(-8, 0, 'PST')
pdt_zone = FixedZone(-7, 0, 'PDT')
pt_zone  = USTimeZone(-8, 0, 'PT', 'PST', 'PDT')

zone_code_map = {
	'UTC': utc_zone,
	'EST': est_zone, 'EDT': edt_zone, 'ET': et_zone,
	'CST': cst_zone, 'CDT': cdt_zone, 'CT': ct_zone,
	'MST': mst_zone, 'MDT': mdt_zone, 'MT': mt_zone,
	'PST': pst_zone, 'PDT': pdt_zone, 'PT': pt_zone}
# - - -  p a r s e A n g l e

def parse_angle(s):
	'''Validate and convert an external angle.

	 [ s is a string ->
		 if s is a valid external angle ->
			return s in radians
		 else -> raise SyntaxError ]
	'''
	#-- 1 --
	minute = second = 0.0
	#-- 2 --
	# [ if s starts with a float followed by 'd' or 'D' ->
	#	 degree	:= that float as type float
	#	 min_tail := s after that float and suffix
	#  else -> raise SyntaxError ]
	degree, min_tail = parse_float_suffix(s, D_PAT, "Degrees followed by 'd'")

	#-- 3 --
	# [ if min_tail is empty -> I
	#  else if min_tail has the form '(float)m' ->
	#	 minute := that (float)
	#  else if min_tail has the form '(float)m(float)s' ->
	#	 minute := the first (float)
	#	 second := the second (float)
	#  else -> raise SyntaxError ]
	if len(min_tail) != 0:
		#-- 3.1 --
		# [ if min_tail starts with a float followed by 'm' or 'M' ->
		#	 minute := that float as type float
		#	 sec_tail := min_tail after all that
		#  else -> raise SyntaxError ]
		minute, sec_tail = parse_float_suffix(min_tail, M_PAT, "Minutes followed by 'm'")

		#-- 3.2 --
		# [ if sec_tail is empty -> I
		#  else if sec_tail starts with a float followed by
		#  's' or 'S' ->
		#	 second := that float as type float
		#	 check_tail := sec_tail after all that
		#  else -> raise SyntaxError ]
		if len(sec_tail) != 0:
			second, check_tail = parse_float_suffix(sec_tail, S_PAT, "Seconds followed by 's'")
			if len(check_tail) != 0:
				raise SyntaxError, ("Unidentifiable angle parts: '{check_tail}'".format(check_tail=check_tail))
	#-- 4 --
	# [ return the angle (degree, minute, second) in radians ]
	angle_degrees = dms_units.mix_to_single((degree, minute, second))
	return radians(angle_degrees)
# - - -  p a r s e F l o a t S u f f i x

def parse_float_suffix(s, code_re, message):
	'''Parse a float followed by a letter code.

	 [ (s is a string) and
		(code_re is a compiled regular expression) and
		(message is a string describing what is expected) ->
		 if s starts with a float, followed by code (using 
		 case-insensitive comparison) ->
			return (x, tail) where x is that float as type float
			and tail is the part of s after the float and code
		 else -> raise SyntaxError, 'Expecting (message)' ]
	'''
	#-- 1 --
	# [ if s starts with a float ->
	#	 x := that float as type float
	#	 code_tail := the part of s after that float
	#  else -> raise SyntaxError, 'Expecting (message)' ]
	x, code_tail = parse_float(s, message)

	#-- 2 --
	# [ if code_tail starts with code (case-insensitive) ->
	#	 return (x, the part of code_tail after the match)
	#  else -> raise SyntaxError ]
	discard, tail = parse_re(code_tail, code_re, message)

	#-- 3 --
	return (x, tail)
# - - -  p a r s e F l o a t

def parse_float(s, message):
	'''Parse a floating-point number at the front of s.

	 [ (s is a string) and
		(message is a string describing what is expected) ->
		 if s begins with a floating-point number ->
			return (x, tail) where x is the number as type float
			and tail is the part of s after the match
		 else -> raise SyntaxError, 'Expecting (message)' ]
	'''
	#-- 1 --
	# [ if the front of s matches FLOAT_PAT ->
	#	 m := a Match object describing the match
	#  else -> raise SyntaxError ]
	raw_float, tail = parse_re(s, FLOAT_PAT, message)

	#-- 2 --
	return (float(raw_float), tail)
# - - -  p a r s e R e

def parse_re(s, regex, message):
	'''Parse a regular expression at the head of a string.

	 [ (s is a string) and
		(regex is a compiled regular expression) and
		(message is a string describing what is expected) ->
		 if s starts with a string that matches regex ->
			return (head, tail) where head is the part of s
			that matched and tail is the rest
		 else ->
			raise SyntaxError, 'Expecting (message)' ]
	'''

	#-- 1 --
	# [ if the head of s matches regex ->
	#	 m := a match object describing the matching part
	#  else -> raise SyntaxError, 'Expecting (message)' ]
	m = regex.match(s)
	if m is None:
		raise SyntaxError, "Expecting {msg}: '{s}'".format(msg=message, s=s)

	#-- 2 --
	# [ return (matched text from s, text from s after match) ]
	head = m.group()
	tail = s[m.end():]
	return (head, tail)
# - - -  p a r s e L a t

def parse_lat(s):
	'''Validate and convert an external latitude.

	 [ s is a nonempty string ->
		 if s is a valid external latitude ->
			return that latitude in radians
		 else -> raise SyntaxError ]
	'''
	#-- 1 --
	# [ last := last character of s
	#  raw_angle := s up to the last character ]
	last = s[-1]
	raw_angle = s[:-1]

	#-- 2 --
	# [ if last matches NS_PAT ->
	#	 ns_flag := last, lowercased
	#  else -> raise SyntaxError ]
	m = NS_PAT.match(last)
	if m is None:
		raise SyntaxError, ("Latitude '{lat}' does not end with 'n' or 's'.".format(lat=s))
	else:
		ns_flag = last.lower()
	#-- 3 --
	# [ if raw_angle is a valid angle ->
	#	 abs_angle := that angle in radians
	#  else -> raise SyntaxError ]
	abs_angle = parse_angle(raw_angle)
	#-- 4 --
	if ns_flag == 's':
		angle = -abs_angle
	else:
		angle = abs_angle

	#-- 5 --
	return angle
# - - -  p a r s e L o n

def parse_lon(s):
	'''Validate and convert an external longitude.

	 [ s is a nonempty string ->
		 if s is a valid external longitude ->
			return that longitude in radians
		 else -> raise SyntaxError ]
	'''
	#-- 1 --
	# [ last := last character of s
	#  raw_angle := s up to the last character ]
	last = s[-1]
	raw_angle = s[:-1]

	#-- 2 --
	# [ if EW_PAT matches last ->
	#	 ew_flag := last, lowercased
	#  else -> raise SyntaxError ]
	m = EW_PAT.match(last)
	if m is None:
		raise SyntaxError, ("Longitude '{lon}' does not end with 'e' or 'w'.".format(lon=s)
	else:
		ew_flag = last.lower()
	#-- 3 --
	# [ if raw_angle is a valid angle ->
	#	 abs_angle := that angle in radians
	#  else -> raise SyntaxError ]
	abs_angle = parse_angle(raw_angle)
	#-- 4 --
	if ew_flag == 'w':
		angle = TWO_PI - abs_angle
	else:
		angle = abs_angle

	#-- 5 --
	return angle
# - - -  p a r s e H o u r s

def parse_hours(s):
	'''Validate and convert a quantity in hours.

	 [ s is a non-empty string ->
		 if s is a valid mixed hours expression ->
			return the value of s as decimal hours
		 else -> raise SyntaxError ]
	'''
	#-- 1 --
	minute = second = 0.0

	#-- 2 --
	# [ if s starts with a float followed by 'h' or 'H' ->
	#	 hour := that float as type float
	#	 min_tail := s after that float and suffix
	#  else -> raise SyntaxError ]
	hour, min_tail = parse_float_suffix(s, H_PAT, "Hours followed by 'h'")

	#-- 3 --
	# [ if min_tail is empty -> I
	#  else if min_tail has the form '(float)m' ->
	#	 minute := that (float)
	#  else if min_tail has the form '(float)m(float)s' ->
	#	 minute := the first (float)
	#	 second := the second (float)
	#  else -> raise SyntaxError ]
	if len(min_tail) != 0:
		#-- 3.1 --
		# [ if min_tail starts with a float followed by 'm' or 'M' ->
		#	 minute := that float as type float
		#	 sec_tail := min_tail after all that
		#  else -> raise SyntaxError ]
		minute, sec_tail = parse_float_suffix(min_tail, M_PAT, "Minutes followed by 'm'")

		#-- 3.2 --
		# [ if sec_tail is empty -> I
		#  else if sec_tail starts with a float followed by
		#  's' or 'S' ->
		#	 second := that float as type float
		#	 check_tail := sec_tail after all that
		#  else -> raise SyntaxError ]
		if len(sec_tail) != 0:
			second, check_tail = parse_float_suffix(sec_tail, S_PAT, "Seconds followed by 's'")
			if len(check_tail) != 0:
				raise SyntaxError, ("Unidentifiable angle parts: '{check_tail}'".format(check_tail=check_tail))
	#-- 4 --
	# [ return the quantity (hour, minute, second) in hours ]
	result = dms_units.mix_to_single((hour, minute, second))
	return result

# - - - - -  c l a s s  M i x e d U n i t s

class MixedUnits(object):
	'''Represents a system with mixed units, e.g., hours/minutes/seconds
	'''
# - - -  M i x e d U n i t s . _ _ i n i t _ _

	def __init__(self, factors):
		'''Constructor
		'''
		self.factors = factors
# - - -  M i x e d U n i t s . m i x T o S i n g l e

	def mix_to_single(self, coeffs):
		'''Convert mixed units to a single value.

		 [ coeffs is a sequence of numbers not longer than
			len(self.factors)+1 ->
			 return the equivalent single value in self's system ]
		'''
		#-- 1 --
		total = 0.0

		#-- 2 --
		# [ if len(coeffs) <= len(self.factors)+1 ->
		#	 coeff_list := a copy of coeffs, right-padded to length
		#		 len(self.factors)+1 with zeroes if necessary ]
		coeff_list = self.__pad(coeffs)
		#-- 3 --
		# [ total +:= (coeff_list[-1] * 
		#		(product of all elements of self.factors)) +
		#	  (coeff_list[-2] *
		#		(product of all elements of self.factors[:-1])) +
		#	  (coeff_list[-3] *
		#		(product of all elements of self.factors[:-2]))
		#		... ]
		for i in range(-1, -len(self.factors)-1, -1):
			total += coeff_list[i]
			total /= self.factors[i]
		#-- 4 --
		total += coeff_list[0]

		#-- 5 --
		return total
# - - -  M i x e d U n i t s . _ _ p a d

	def __pad(self, coeffs):
		'''Pad coefficient lists to standard length.

		 [ coeffs is a sequence of numbers ->
			 if len(coeffs) > len(self.factors)+1 ->
				raise ValueError
			 else ->
				return a list containing the elements of coeff,
				plus additional zeroes on the right if necessary
				so that the result has length len(self.factors)+1 ]
		'''
		#-- 1 --
		# [ std_len := 1 + len(self.factors)
		#  shortage := 1 + len(self.factors) - len(coeffs)
		#  result := a copy of coeffs as a list ]
		std_len = 1 + len(self.factors)
		shortage = std_len - len(coeffs)
		result = list(coeffs)

		#-- 2 --
		# [ if shortage < 0 ->
		#	 raise ValueError
		#  else ->
		#	 result := result + (a list of shortage zeroes) ]
		if shortage < 0:
			raise ValueError, ('Value {coeffs} has too many elements; max is {std_len}.'.format(coeffs=coeffs, std_len=std_len))
		elif shortage > 0:
			result += [0.0] * shortage

		#-- 3 --
		return result
# - - -  M i x e d U n i t s . s i n g l e T o M i x

	def single_to_mix(self, value):
		'''Convert to mixed units.

		 [ value is a float ->
			 return value as a sequence of coefficients in
			 self's system ]
		'''
		#-- 1 --
		# [ whole := whole part of value
		#  frac := fractional part of value ]
		whole, frac = divmod(value, 1.0)
		result = [int(whole)]

		#-- 2 --
		# [ result := result with integral parts of value
		#			  in self's system appended ]
		for factorx in range(len(self.factors)):
			frac *= self.factors[factorx]
			whole, frac = divmod(frac, 1.0)
			result.append(int(whole))

		#-- 3 --
		# [ result := result with frac added to its last element ]
		result[-1] += frac

		#-- 4 --
		return result
# - - -  M i x e d U n i t s . f o r m a t

	def format(self, coeffs, decimals=0, lz=False):
		'''Format mixed units.

		 [ (coeffs is a sequence of numbers as returned by
			MixedUnits.single_to_mix()) and
			(decimals is a nonnegative integer) and
			(lz is a bool) ->
			 return a list of strings corresponding to the values
			 of coeffs, with all the values but the last formatted
			 as integers, all values zero padded iff lz is true,
			 and the last value with (decimals) digits after the
			 decimal point ]
		'''
		#-- 1 --
		coeff_list = self.__pad(coeffs)

		#-- 2 --
		# [ result := the values from coeff_list[:-1] formatted
		#			  as integers ]
		if lz:
			fmt = '%02d'
		else:
			fmt = '%d'
		result = [fmt % x for x in coeff_list[:-1]]

		#-- 2 --
		# [ whole := whole part of coeff_list[-1]
		#  frac  := fractional part of coeff_list[-1]
		#  fuzz  := 0.5 * (10 ** (-decimals) ]
		whole, frac = divmod(float(coeff_list[-1]), 1.0)
		fuzz = 0.5 * (10.0 ** (-decimals))

		#-- 3 --
		# [ if frac >= (1-fuzz) ->
		#	 result +:= [whole+frac-fuzz], formatted with
		#				 (decimals) digits after the decimal
		#  else ->
		#	 result +=  coeff_list[-1], formatted with (decimals)
		#				 digits after the decimal ]
		if frac >= (1.0 - fuzz):
			corrected = whole + frac - fuzz
		else:
			corrected = coeff_list[-1]

		#-- 4 --
		# [ if lz ->
		#	 s := corrected, formatted with 2 digits of left-zero
		#			padding and (decimals) precision
		#  else ->
		#	 s := corrected, formatted with (decimals) precision ]
		if lz:
			if decimals:
				n = decimals + 3
			else:
				n = decimals + 2

			s = '%0*.*f' % (n, decimals, corrected)
		else:
			s = '%.*f' % (decimals, corrected)
			 
		#-- 5 --
		result.append(s)

		#-- 6 --
		return result
dms_units = MixedUnits((60, 60))

# - - - - -  c l a s s  L a t L o n

class LatLon(object):
	'''Represents a latitude+longitude.
	'''
# - - - L a t L o n . _ _ i n i t _ _

	def __init__(self, lat, lon):
		'''Constructor for LatLon.
		'''
		self.lat = lat
		self.lon = lon % TWO_PI
# - - -  L a t L o n . _ _ s t r _ _

	def __str__(self):
		'''Return self as a string.
		'''
		#-- 1 --
		if self.lon >= math.pi:
			e_w = 'W'
			lon_deg = degrees(TWO_PI - self.lon)
		else:
			e_w = 'E'
			lon_deg = degrees(self.lon)

		#-- 2 --
		if self.lat < 0:
			n_s = 'S'
			lat_deg = degrees(-self.lat)
		else:
			n_s = 'N'
			lat_deg = degrees(self.lat)
		#-- 3 --
		# [ lat_list := three formatted values of lat_deg in
		#				degrees/minutes/seconds
		#  lon_list := three formatted values of lon_deg similarly ]
		lat_list = dms_units.format(dms_units.single_to_mix(lat_deg), 1)
		lon_list = dms_units.format(dms_units.single_to_mix(lon_deg), 1)

		#-- 4 --
		return '[{lat_d}d {lat_s}\' {lat_m}" {ns} Lat {lon_d}d {lon_s}\' {lon_m}" {ew} Lon]'.format(lat_d=lat_list[0], lat_s=lat_list[1],
																									lat_m=lat_list[2], ns=n_s,
																									lon_d=lon_list[0], lon_s=lon_list[1],
																									lon_m=lon_list[2], ew=e_w)

# - - - - -  c l a s s  J u l i a n D a t e

class JulianDate(object):
	'''Class to represent Julian-date timestamps.

	 State/Invariants:
		.f: [ (Julian date as a float) - JULIAN_BIAS ]
	'''
# - - -  J u l i a n D a t e . _ _ i n i t _ _

	def __init__(self, j, f=0.0):
		'''Constructor for JulianDate.
		'''
		self.j = j - JULIAN_BIAS + f
# - - -  J u l i a n D a t e . _ _ f l o a t _ _

	def __float__(self):
		'''Convert self to a float.
		'''
		return self.j + JULIAN_BIAS
# - - -  J u l i a n D a t e . d a t e t i m e

	def datetime(self):
		'''Convert to a standard Python datetime object in UT.
		'''
		#-- 1 --
		# [ i := int(self.j + 0.5)
		#  f := (self.j + 0.5) % 1.0 ]
		i, f = divmod(self.j + 0.5, 1.0)
		i += JULIAN_BIAS
		#-- 2 --
		if i > 2299160:
			a = int((i - 1867216.25) / 36524.25)
			b = i + 1 + a - int(a / 4.0)
		else:
			b = i
		#-- 3 --
		c = b + 1524
		#-- 4 --
		d = int((c - 122.1) / 365.25)
		#-- 5 --
		e = int(365.25 * d)
		#-- 6 --
		g = int((c - e) / 30.6001)
		#-- 7 --
		day_frac = c - e + f - int(30.6001 * g)
		day, frac = divmod(day_frac, 1.0)
		dd = int(day)
		hr, mn, sc = dms_units.single_to_mix(24.0 * frac)
		#-- 8 --
		if g < 13.5:
			mm = int(g - 1)
		else:
			mm = int(g - 13)
		#-- 9 --
		if mm > 2.5:
			yyyy = int(d - 4716)
		else:
			yyyy = int(d - 4715)
		#-- 10 --
		sec, frac_sec = divmod(sc, 1.0)
		usec = int(frac_sec * 1e6)
		return datetime.datetime(yyyy, mm, dd, hr, mn, int(sec), usec)

# - - -  J u l i a n D a t e . o f f s e t

	def offset(self, delta):
		'''Return a new JulianDate for self+(delta days)

		 [ delta is a number of days as a float ->
			 return a new JulianDate (delta) days in the
			 future, or past if negative ]
		'''
		#-- 1 --
		new_j = self.j + delta

		#-- 2 --
		# [ new_whole := whole part of new_j
		#  new_frac  := fractional part of new_j ]
		new_whole, new_frac = divmod(new_j)

		#-- 3 --
		return JulianDate(new_whole + JULIAN_BIAS, new_frac)
# - - -  J u l i a n D a t e . _ _ s u b _ _

	def __sub__(self, other):
		'''Implement subtraction.

		 [ other is a JulianDate instance ->
			 return self.j - other.j ]
		'''
		return self.j - other.j
# - - -  J u l i a n D a t e . _ _ c m p _ _

	def __cmp__(self, other):
		'''Compare two instances.

		 [ other is a JulianDate instance ->
			 if self.j < other.j -> return a negative number
			 else if self.j == other.j -> return zero
			 else -> return a positive number ]
		'''
		return cmp(self.j, other.j)
# - - -  J u l i a n D a t e . f r o m D a t e t i m e

#  @staticmethod
	def from_datetime(dt):
		'''Create a JulianDate instance from a datetime.datetime.

		 [ dt is a datetime.datetime instance ->
			 if dt is naive ->
				return the equivalent new JulianDate instance,
				assuming dt expresses UTC
			 else ->
				return a new JulianDate instance for the UTC
				time equivalent to dt ]			 
		'''
		#-- 1 --
		# [ if dt is naive ->
		#	 utc := dt
		#  else ->
		#	 utc := dt - dt.utc_offset() ]
		utc = dt
		offset = dt.utc_offset()
		if offset:
			utc = dt - offset
		#-- 2 --
		# [ frac_day := fraction of a day in [0.0,1.0) made from
		#	  utc.hour, utc.minute, utc.second, and utc.microsecond ]
		s = float(utc.second) + float(utc.microsecond) * 1e-6
		hours = dms_units.mix_to_single((utc.hour, utc.minute, s))
		frac_day = hours / 24.0
		#-- 3 --
		y = utc.year
		m = utc.month
		d = utc.day
		#-- 4 --
		if m <= 2:
			y, m = y-1, m+12
		#-- 5 --
		if (y, m, d) >= (1582, 10, 15):
			A = int(y / 100)
			B = 2 - A + int(A / 4)
		else:
			B = 0
		#-- 6 --
		C = int(365.25 * y)
		D = int(30.6001 * (m + 1))
		#-- 7 --
		# [ if frac_day+0.5 >= 1.0 ->
		#	 s += 1
		#	 frac_day := (frac_day + 0.5) % 1.0
		#  else ->
		#	 frac_day := frac_day + 0.5 ]
		day_carry, frac_day = divmod(frac_day + 0.5, 1.0)
		d += day_carry

		#-- 8 --
		j = B + C + D + d + 1720994

		#-- 9 --
		return JulianDate(j, frac_day)
	from_datetime = staticmethod(from_datetime)

# - - - - -  c l a s s  S i d e r e a l T i m e

class SiderealTime(object):
	'''Represents a sidereal time value.

	 State/Internals:
		.hours:	 [ self as 15-degree hours ]
		.radians:  [ self as radians ]
	'''
# - - -  S i d e r e a l T i m e . _ _ i n i t _ _

	def __init__(self, hours):
		'''Constructor for SiderealTime
		'''
		self.hours = hours % 24.0
		self.radians = hours_to_radians(self.hours)
# - - -  S i d e r e a l T i m e . _ _ s t r _ _

	def __str__(self):
		'''Convert to a string such as '[04h 40m 5.170s]'.
		'''

		#-- 1 --
		# [ values := self.hours as a list of mixed units
		#	  in dms_units terms, formatted as left-zero
		#	  filled strings with 3 digits after the decimal ]
		mix = dms_units.single_to_mix(self.hours)
		values = dms_units.format(mix, decimals=3, lz=True)

		#-- 2 --
		values = tuple(values)
		return '[{hh}h {mm}m {ss}s]'.format(hh=values[0], mm=values[1], ss=values[2])
# - - -  S i d e r e a l T i m e . u t c

	def utc(self, date):
		'''Convert GST to UTC.

		 [ date is a UTC date as a datetime.date instance ->
			 return the first or only time at which self's GST
			 occurs at longitude 0 ]
		'''
		#-- 1 --
		# [ n_days := number of days between Jan. 0 of year
		#	  (date.year) and date ]
		n_days = day_num(date)

		#-- 2 --
		# [ t0 := (n_days * A - B(date.year)), normalized to
		#		  interval [0,24) ]
		t0 = (n_days * SIDEREAL_A - SiderealTime.factor_b(date.year)) % 24.0

		#-- 3 --
		# [ t1 := ((self in decimal hours)-t0), normalized to
		#		  the interval [0,24) ]
		t1 = (radians_to_hours(self.radians) - t0) % 24.0

		#-- 4 --
		gmt_hours = t1 * 0.997270

		#-- 5 --
		# [ dt := a datetime.datetime instance whose date comes
		#		  from (date) and whose time is (gmt_hours)
		#		  decimal hours ]
		hour, minute, float_sec = dms_units.single_to_mix(gmt_hours)
		whole_sec, frac_sec = divmod(float_sec, 1.0)
		second = int(whole_sec)
		micros = int(frac_sec * 1e6)
		dt = datetime.datetime(date.year, date.month, date.day, hour, minute, second, micros)
		
		#-- 6 --
		return dt
# - - -  S i d e r e a l T i m e . f a c t o r B

#  @staticmethod
	def factor_b(yyyy):
		'''Compute sidereal conversion factor B for a given year.

		 [ yyyy is a year number as an int ->
			 return the GST at time yyyy-01-00T00:00 ]
		'''
		#-- 1 --
		# [ jan_jd := the Julian date of January 0.0 of year
		#			 (yyyy), as a float ]
		jan_dt = datetime.datetime(yyyy, 1, 1)
		jan_jd = float(JulianDate.from_datetime(jan_dt)) - 1.0

		#-- 2 --
		s = jan_jd - 2415020.0

		#-- 3 --
		t = s / 36525.0

		#-- 4 --
		r = (0.00002581 * t + 2400.051262) * t + 6.6460656

		#-- 5 --
		u = r - 24 * (yyyy - 1900)

		#-- 6 --
		return 24.0 - u

	factor_b = staticmethod(factor_b)
# - - -  S i d e r e a l T i m e . g s t

	def gst(self, e_lon):
		'''Convert LST to GST.

		 [ self is local sidereal time at longitude e_lon
			radians east of Greenwich ->
			 return the equivalent GST as a SiderealTime instance ]
		'''
		#-- 1 --
		# [ delta_hours := e_lon expressed in hours ]
		delta_hours = radians_to_hours(e_lon)

		#-- 2 --
		gst_hours = (self.hours - delta_hours) % 24.0

		#-- 3 --
		return SiderealTime(gst_hours)
# - - -  S i d e r e a l T i m e . l s t

	def lst(self, e_lon):
		'''Convert GST to LST.

		 [ (self is Greenwich sidereal time) and
			(e_lon is a longitude east of Greenwich in radians) ->
			 return a new SiderealTime representing the LST
			 at longitude e_lon ]
		'''
		#-- 1 --
		# [ delta_hours := e_lon expressed in hours ]
		delta_hours = radians_to_hours(e_lon)

		#-- 2 --
		gmt_hours = (self.hours + delta_hours) % 24.0

		#-- 3 --
		return SiderealTime(gmt_hours)
# - - -  S i d e r e a l T i m e . f r o m D a t e t i m e

	SIDEREAL_C = 1.002738

#  @staticmethod
	def from_datetime(dt):
		'''Convert civil time to Greenwich Sidereal.

		 [ dt is a datetime.datetime instance ->
			 if dt has time zone information ->
				return the GST at the UTC equivalent to dt
			 else ->
				return the GST assuming dt is UTC ]
		'''
		#-- 1 --
		# [ if dt is naive ->
		#	 utc := dt
		#  else ->
		#	 utc := the UTC time equivalent to dt ]
		utc = dt
		tz = dt.tzinfo
		if tz is not None:
			offset = tz.utc_offset(dt)
			if offset is not None:
				utc = dt - offset

		#-- 2 --
		# [ n_days := number of days between January 0.0 and utc ]
		n_days = day_num(utc)

		#-- 3 --
		t0 = n_days * SIDEREAL_A - SiderealTime.factor_b(utc.year)
		#-- 4 --
		# [ dec_utc := utc as decimal hours ]
		float_sec = utc.second + float(utc.microsecond) / 1e6
		dec_utc = dms_units.mix_to_single((utc.hour, utc.minute, float_sec))

		#-- 5 --
		# [ gst := (dec_utc * C + t0), normalized to interval [0,24) ]
		gst = (dec_utc * SiderealTime.SIDEREAL_C + t0) % 24.0

		#-- 6 --
		return SiderealTime(gst)
	from_datetime = staticmethod(from_datetime)

# - - - - -  c l a s s  A l t A z

class AltAz(object):
	'''Represents a sky location in horizon coords. (altitude/azimuth)

	 Exports/Invariants:
		.alt:  [ altitude in radians, in [-math.pi, +math.pi] ]
		.az:	[ azimuth in radians, in [0, 2 * math.pi] ]
	'''
# - - -  A l t A z . _ _ i n i t _ _

	def __init__(self, alt, az):
		'''Constructor for AltAz, horizon coordinates.

		 [ (alt is an altitude in radians) and
			(az is an azimuth in radians) ->
			 return a new AltAz instance with those values,
			 normalized as per class invariants ]
		'''
		self.alt = alt
		self.az = az
# - - -  A l t A z . r a D e c

	def ra_dec(self, lst, lat_lon):
		'''Convert horizon coordinates to equatorial.

		 [ (lst is a local sidereal time as a SiderealTime instance) and
			(lat_lon is the observer's position as a LatLon instance) ->
			 return the corresponding equatorial coordinates as a
			 RADec instance ]			
		'''
		#-- 1 --
		# [ dec := declination of self at lat_lon in radians
		#  hour_rads := hour angle of self at latlon in radians ]
		dec, hour_rads = coord_rotate(self.alt, lat_lon.lat, self.az)

		#-- 2 --
		# [ hour_rads is an hour angle in radians ->
		#	 h := hour_rads in hours ]
		h = radians_to_hours(hour_rads)
		#-- 3 --
		# [ ra := right ascension for hour angle (h) at local
		#		  sidereal time (lst) and location (lat_lon) ]
		ra = hours_to_radians((lst.hours - h) % 24.0)

		#-- 4 --
		return RADec(ra, dec)
# - - -  A l t A z . _ _ s t r _ _

	def __str__(self):
		'''Convert self to a string.
		'''
		#-- 1 --
		# [ alt_list := self.alt, formatted as degrees, minutes,
		#	  and seconds
		#  az_list := self.az, formatted as degrees, minutes, and
		#	  seconds ]
		alt_list = dms_units.format(dms_units.single_to_mix(degrees(self.alt)), lz=True, decimals=3)
		az_list = dms_units.format(dms_units.single_to_mix(degrees(self.az)), lz=True, decimals=3)

		#-- 2 --
		az_list = tuple(az_list)
		alt_list = tuple(alt_list)
		return '[az {az_d}d {az_m}\' {az_s}" alt {alt_d}d {alt_m}\' {alt_s}"]'.format(az_d=az_list[0], az_m=az_list[1], az_s=az_list[2],
																						alt_d=alt_list[0], alt_m=alt_list[1], alt_s=alt_list[2])
# - - -  c o o r d R o t a t e

def coord_rotate(x, y, z):
	'''Used to convert between equatorial and horizon coordinates.

	 [ x, y, and z are angles in radians ->
		 return (xt, yt) where
		 xt = arcsin(sin(x) * sin(y) + cos(x) * cos(y) * cos(z)) and
		 yt = arccos((sin(x) - sin(y) * sin(xt)) / (cos(y) * cos(xt))) ]
	'''
	#-- 1 --
	xt = asin(sin(x) * sin(y) + cos(x) * cos(y) * cos(z))

	#-- 2 --
	yt = acos((sin(x) - sin(y) * sin(xt)) / (cos(y) * cos(xt)))

	#-- 3 --
	if sin(z) > 0.0:
		yt = TWO_PI - yt

	#-- 4 --
	return (xt, yt)

# - - - - -  c l a s s  R A D e c

class RADec(object):
	'''Represents a celestial location in equatorial coordinates.

	 Exports/Invariants:
		.ra:	 [ right ascension in radians ]
		.dec:	 [ declination in radians ]
	'''
# - - -  R A D e c . _ _ i n i t _ _

	def __init__(self, ra, dec):
		'''Constructor for RADec.
		'''
		self.ra = ra % TWO_PI
		self.dec = dec
# - - -  R A D e c . h o u r A n g l e

	def hour_angle(self, utc, e_lon):
		'''Find the hour angle at a given observer's location.

		 [ (utc is a Universal Time as a datetime.datetime) and
			(e_lon is an east longitude in radians) ->
			 return the hour angle of self at that time and
			 longitude, in radians ]
		'''
		return ra_to_hour_angle(self.ra, utc, e_lon)
# - - -  R A D e c . a l t A z

	def alt_az(self, h, lat):
		'''Convert equatorial to horizon coordinates.

		 [ (h is an object's hour angle in radians) and
			(lat is the observer's latitude in radians) ->
			 return self's position in the observer's sky
			 in horizon coordinates as an AltAz instance ]
		'''
		#-- 1 --
		# [ alt := altitude of self as seen from lat_lon at utc
		#  az := azimuth of self as seen from lat_lon at utc ]
		alt, az = coord_rotate(self.dec, lat, h)

		#-- 2 --
		return AltAz(alt, az)
# - - -  R A D e c . _ _ s t r _ _

	def __str__(self):
		'''Return self as a string.
		'''
		#-- 1 --
		# [ ra_units := units of self.ra as hours/minutes/seconds
		#  dec_units := units of self.dec as degrees/minutes/seconds
		ra_units = dms_units.format(dms_units.single_to_mix(radians_to_hours(self.ra)), lz=True, decimals=3)
		dec_units = dms_units.format(dms_units.single_to_mix(degrees(self.dec)), lz=True, decimals=3)

		#-- 2 --
		ra_units = tuple(ra_units)
		dec_units = tuple(dec_units)
		return '[{ra_h}h {ra_m}m {ra_s}s, {dec_d}d {dec_m}\' {dec_s}"]'.format(ra_d=ra_units[0], ra_m=ra_units[1], ra_s=ra_units[2],
																				dec_d=dec_units[0], dec_m=dec_units[1], dec_s=dec_units[2])
