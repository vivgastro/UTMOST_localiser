#!/home/observer/miniconda2/bin/python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pablo_coords as pc
import subprocess as sp
import re, sys, os, argparse, ephem
from datetime import datetime as DT
from datetime import timedelta as TD
from parse_cfg import parse_cfg as pcfg
from astropy.coordinates import SkyCoord, FK5, CIRS
from astropy.time import Time
import astropy.units as u
#import scipy.interpolate, scipy.optimize
from scipy.optimize import curve_fit


#test_psr_coords = ['19:35:47.83', '+16:16:39.98']
#test_psr_coords = ['16:44:49.281', '-45:59:09.5']
test_psr_coords = ['12:43:17.158', '-64:23:23.85']
#test_psr_coords = ['18:14:41.26', '-06:18:01.7']
#test_psr_coords = ['19:32:14.06', '+10:59:33.38']
#test_psr_coords = ['05:34:31.973', '+22:00:52.06']

utc_start_format = "%Y-%m-%d-%H:%M:%S"
utc_format = "%Y-%m-%d-%H:%M:%S.%f"

fmt_str = lambda(utc_str): utc_format if '.' in utc_str else utc_start_format
fmt_obj = lambda(utc_obj): utc_format if '.' in utc_obj.isoformat() else utc_start_format

def radec_to_J2000(ra, dec, utc):
  eq = ephem.Equatorial(ra, dec, epoch = utc)
  eq2 = ephem.Equatorial(eq, epoch=ephem.J2000)
  return eq2

def radecJ2000_to_utc(ra, dec, utc):
  coords = "{0} {1}".format(ra, dec)
  ap_c = SkyCoord(coords, unit=(u.hourangle, u.deg), frame=FK5(equinox="J2000"))
  x = DT(utc, fmt_str(utc))
  epoch = Time(x)
  new_coords = ap_c.transform(FK5(equinox=epoch))

def get_lst(utc):
  x = utc.strftime(fmt_obj(utc)).split(".")
  if len(x)==1:
    cmd = "mopsr_getlst {0}".format(x[0])
  elif len(x)==2:
    cmd = "mopsr_getlst {0} -f 0.{1}".format(x[0], x[1])
  else:
    raise ValueError("Invalid UTC provided {}".format(utc.strftime(fmt_obj(utc))))
  stdout, stderr = sp.Popen(cmd, shell=True, stdout=sp.PIPE).communicate()
  patt = r'radians\=\d\.\d*'
  x = re.compile(patt)
  val = float(x.findall(stdout)[0].strip('radians='))
  return np.rad2deg(val)

def get_ut1(utc_start):
  path = "/data/mopsr/results/{0}/FB/obs.header".format(utc_start)
  if not os.path.exists(path):
    path = "/data/mopsr/old_results/{0}/FB/obs.header".format(utc_start)
    if not os.path.exists(path):
      raise IOError("File not found: {0} and UT1 offset not provided".format(path))

  obsheader = pcfg(path)
  ut1 = obsheader.UT1_OFFSET
  return ut1

def correct_for_ut1(utc, ut1):
  if ut1==None:
    ut1 = get_ut1(utc.strftime(fmt_obj(utc)))
  ut1 = TD(seconds=ut1)
  print ut1.total_seconds()
  cutc = utc + ut1
  return cutc

def get_utc_pulse(utc_start, toff):
  toff = TD(seconds = toff)
  return utc_start + toff

def get_MD(pulse):
  center_fb = int((pulse['nbeams']+1)/2.0) + 1
  md = (pulse['fb'] - center_fb) * pulse['fbs']
  md /= np.cos(np.deg2rad(md))      #Check this, can be wrong or should be removed
  return md

def get_lst_pulse(pulse):
  utc_start = DT.strptime(pulse['utc_start'], fmt_str(pulse['utc_start']))
  utc_start_corrected = correct_for_ut1(utc_start, args.ut1)
  utc_pulse = get_utc_pulse(utc_start_corrected, pulse['tstamp'])
  lst_pulse = get_lst(utc_pulse)
  print "UTC_pulse = {0}".format(utc_pulse)
  #pulse['utc_pulse'] = utc_pulse.strftime(fmt_obj(utc_pulse))
  print "LST = {0}".format(lst_pulse)
  
  return lst_pulse

#@profile
def get_RA_DEC_values(NS, MD, pulse):
  HA = np.rad2deg(pc.HA(NS, MD))
  DEC = np.rad2deg(pc.DEC(NS, MD))

  RA_zenith = get_lst_pulse(pulse)
  RA = RA_zenith - HA

  return RA, DEC

def to_ICRS(RA, DEC, utc):
  utc = DT.strptime(utc, fmt_str(utc))
  coords = SkyCoords(RA, DEC, unit=(u.deg, u.deg), frame=FK5(equinox=utc))
  coords_ICRS = coords.transform_to(frame=ICRS)
  return coords_ICRS.ra.value, coords_ICRS.dec.value

#@profile
def get_RA_DEC_limits(NS, MD, pulse):
  lower_MD = MD - np.deg2rad(pulse['fbs']/4)
  upper_MD = MD + np.deg2rad(pulse['fbs']/4)
  RA, DEC = get_RA_DEC_values(NS, MD, pulse)
  RA = RA[1:-1]
  DEC = DEC[1:-1]
  RA_l, DEC_l = get_RA_DEC_values(NS, lower_MD, pulse)
  #RA_lower = scipy.interpolate.interp1d(DEC_l, RA_l, kind='linear')(DEC)
  RA_lower = np.interp(DEC, DEC_l, RA_l)

  RA_u, DEC_u = get_RA_DEC_values(NS, upper_MD, pulse)
  #RA_upper = scipy.interpolate.interp1d(DEC_u, RA_u, kind='cubic')(DEC)
  RA_upper = np.interp(DEC, DEC_u, RA_u)

  return DEC, RA_lower, RA, RA_upper, DEC_l, DEC_u, RA_l, RA_u

def get_NS0_DEC0(pulse):
  a = "/data/mopsr/results/{}/obs.info".format(pulse['utc_start'])
  b = "/data/mopsr/old_results/{}/obs.info".format(pulse['utc_start'])
  obs_info_file = a if os.path.exists(a) else b

  x = pcfg(obs_info_file)
  DEC0 = x.DEC.split(":")
  sign = np.sign(float(DEC0[0]))
  DEC0 = np.abs(float(DEC0[0])) + float(DEC0[1])/60. + float(DEC0[2])/3600.
  DEC0 *= sign

  NS0 = DEC0 - np.rad2deg(pc.lat)
  return NS0, DEC0

def fit_around_DEC0(model, RA, DEC, DEC0):
  fit_params, fit_covars = curve_fit(model, DEC, RA, p0=[0, 0, 0, 0])
  return fit_params

#@profile
def plot(pulse, ax, idx):
  NS0, DEC0 = get_NS0_DEC0(pulse)
  MD = np.deg2rad(get_MD(pulse))
  npoints = args.NS_ex * 1e2  #100 points per degree
  NS = np.deg2rad(np.linspace(NS0 - args.NS_ex, NS0 + args.NS_ex, npoints))
  
  DEC, RA_lower, RA, RA_upper, DEC_l, DEC_u, RA_l, RA_u = get_RA_DEC_limits(NS, MD, pulse)
  cmap = mpl.cm.Set1
  if idx ==1:
    model = lambda dec, a, b, c, d: a + b*(dec-DEC0) + c*(dec-DEC0)**2 + d*(dec-DEC0)**3
    fit = fit_around_DEC0(model, RA, DEC, DEC0)
    ax.plot(model(DEC, *fit)/15, DEC, 'k--' )
  
  ax.fill_betweenx(DEC, RA_lower/15, RA_upper/15, color=cmap(0%cmap.N), alpha = 0.1, label='Tstamp:{0} FB:{1}'.format(pulse['tstamp'], pulse['fb']) )
  ax.plot(RA/15, DEC, '-', c = cmap(idx%cmap.N)  )
  ax.plot(RA_l/15, DEC_l, '--', c = cmap(idx%cmap.N) )
  ax.plot(RA_u/15, DEC_u, '--', c = cmap(idx%cmap.N) )



#@profile
def plot_test_psr(ax, pulse):
  coords = SkyCoord("{0} {1}".format(test_psr_coords[0], test_psr_coords[1]), unit=(u.hourangle, u.deg), frame=FK5(equinox="J2000"))
  coords_today = coords.transform_to(FK5(equinox=Time(DT.strptime(pulse['utc_start'], fmt_str(pulse['utc_start'])))))

  ax.plot(coords.ra.value/15, coords.dec.value, 'r+', label="Pulsar J2000")
  ax.plot(coords_today.ra.value/15, coords_today.dec.value, 'bo', label="Puslar now")

def plot_test_psr_correcting_for_abberation(ax, pulse):
  from astropy._erfa import eo06a
  from astropy.coordinates.builtin_frames.utils import get_jd12

  coords = SkyCoord("{0} {1}".format(test_psr_coords[0], test_psr_coords[1]), unit=(u.hourangle, u.deg), frame=FK5(equinox="J2000"))
  time = Time(DT.strptime(pulse['utc_start'], fmt_str(pulse['utc_start'])))
  cirs = coords.transform_to(CIRS(obstime=time))

  eo = eo06a(*get_jd12(time, 'tt')) * u.rad
  correct_coords = ((cirs.ra - eo).value, cirs.dec.value)
  ax.plot(correct_coords[0]/15., correct_coords[1], 'ro', label="Puslar abberated")

#@profile
def main(args):
  #Read the pulses file in a numpy structured array
  dtype = np.dtype([('utc_start', 'S19'), ('tstamp', np.float), ('fb', np.int16), ('nbeams', np.int16), ('fbs', np.float)])
  pulses = np.loadtxt(args.pfile, dtype=dtype)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  for ii, pulse in enumerate(pulses):
    plot(pulse, ax, ii)
  
  plot_test_psr(ax, pulses[0])
  plot_test_psr_correcting_for_abberation(ax, pulses[0])
  plt.legend()
  plt.xlabel("RA")
  plt.ylabel("DEC")
  plt.show()
  #plt.savefig("/home/vgupta/tmp/tmp.png")

if __name__ == '__main__':
  a = argparse.ArgumentParser()
  a.add_argument("pfile", type=str, help="Path to the file containing info about the  pulses")
  a.add_argument("-ut1", type=float, help="UT1 offset in seconds (def:Will be taken from obs.header)", default=None)
  a.add_argument("-NS_ex", type=float, help="North-South extent to which draw the loc arc on either side (def:20 deg)", default=20)

  args= a.parse_args()
  main(args)
