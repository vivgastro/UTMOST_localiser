from __future__ import print_function
import pablo_coords as pc
from astropy.coordinates import SkyCoord, Angle, FK5, CIRS
from astropy._erfa import eo06a
import astropy.units as u
from astropy.coordinates.builtin_frames.utils import get_jd12
import numpy as np
from astropy.time import Time 
from datetime import datetime as DT
from datetime import timedelta as TD
import argparse, os

from parse_cfg import parse_cfg as pcfg

utc_start_format = "%Y-%m-%d-%H:%M:%S"
utc_format = "%Y-%m-%d-%H:%M:%S.%f"
#mol_long = Angle('149.42465783333333d')   #from mopsr_delays.h, only affects results by ~1 arc-second on sky (0.1 seconds in RA)
mol_long = Angle('149.425d')

fmt_str = lambda(utc_str): utc_format if '.' in utc_str else utc_start_format
fmt_obj = lambda(utc_obj): utc_format if '.' in utc_obj.isoformat() else utc_start_format

OLD_RESULTS_PATH = "/data/mopsr/old_results"
RESULTS_PATH = "/data/mopsr/results"

def read_obs_header(utc_start = None, obs_header_file = None):
  
  if obs_header_file:
    if os.path.exists(obs_header_file):
      obs_header_path = obs_header_file
    else:
      raise RuntimeError("obs.header file - {0} does not exist".format(obs_header_file))
  
  else:
    if utc_start is None:
      raise RuntimeError("Need to specify either utc_start or obs_header_file")

    old_utc_path = os.path.join(OLD_RESULTS_PATH, utc_start)
    new_utc_path = os.path.join(RESULTS_PATH, utc_start)

    if os.path.exists(old_utc_path):
      utc_path = old_utc_path
    elif os.path.exists(new_utc_path):
      utc_path = new_utc_path
    else:
      raise RuntimeError("Could not find {0} in {1} or {2}".format(utc_start, RESULTS_PATH, OLD_RESULTS_PATH))

    obs_header_path = os.path.join(utc_path, "FB/obs.header")

  if not os.path.exists(obs_header_path):
    raise RuntimeError("Could not find {}".format(obs_header_path))
  obs_header = pcfg(obs_header_path)

  return obs_header



def from_J2000_to_apparent(RA, DEC, utc):
  #RA and DEC need to be in degrees. Both can be arrays
  coords = SkyCoord(RA, DEC, frame=FK5(equinox="J2000"))
  time = Time(utc)
  cirs = coords.transform_to(CIRS(obstime=time))

  eo = eo06a(*get_jd12(time, 'tt')) * u.rad
  correct_ra = (cirs.ra - eo)
  correct_dec = cirs.dec
  #correct_coords = ((cirs.ra - eo).value, cirs.dec.value)
  return correct_ra, correct_dec
  #return correct_coords.ra.deg, correct_coords.dec.deg

def get_lst(utc):
  lst = Time(utc).sidereal_time(kind='apparent', longitude=mol_long)
  return lst

def get_HA_rad(RA, DEC, utc):
  lst = Angle(get_lst(utc), unit=u.hourangle)
  #Taking convention -- HA is positive towards east
  HA = RA - lst     #If RA < lst (i.e. the source has transitted), then HA will be negative -- which is technically incorrect
                    #But this HA will be passed to pc.MD and pc.NS, which only use HA inside sine and cosine functions
                    #And sin(-theta) = sin(360-theta), therefore I am happy to keep HA value in negative here

  HA_rad = HA.rad * np.cos(DEC.rad)
  return HA_rad

def get_MD_from_RA_DEC(RA, DEC, utc):
  '''
  Computes the MD of a source at a given time
  Requires the RA, DEC (J2000) of the source, and the utc at which you want to determine the MD location of the source

  Input
  -----
  RA:     RA of the source (astropy.Angle object)
  DEC:    DEC of the source (astropy.Angle object)
  utc:    datetime.datetime object

  Output
  ------
  MD:     MD of the source in radians (astropy.Angle object)
  '''
  
  RA, DEC = from_J2000_to_apparent(RA, DEC, utc)
  HA_rad = get_HA_rad(RA, DEC, utc)

  #NS = Angle(pc.NS(HA_rad, DEC.rad), unit=u.rad)
  #MD = Angle(pc.MD(HA_rad, DEC.rad), unit=u.rad)

  MD0_rad = pc.MD(0, DEC.rad)
  MD_pulse_rad =  Angle(HA_rad - MD0_rad, unit=u.rad)
  
  return  MD_pulse_rad


def get_fanbeam_from_coordinates(RA, DEC, utc_start=None, tcand=None, utc_cand=None, fb_s = None, cfb = None, obs_header_file = None):
  '''
  Computes the fanbeam location of a source at a given utc.
  
  Input
  -----
  RA:     RA of the source (astropy.Angle object -- can be a list or an array) or (str -- only one value)
  DEC:    DEC of the source (astropy.Angle object -- can be a list or an array) or (str -- only one value)

  One of the following combinations:

  utc_start:    UTC_START of the observation (str)
  tcand:        Tcand of the candidate (float)

  OR

  utc_cand:     UTC of the candidate (str)

  One of the following combinations:

  utc_start:    UTC_START of the observation (str)
   OR
  fb_s:         Fan-beam spacing in degrees (float)
                If fb_s is not provided, the function tries to read it from the UTC_START

  One of the following combinatios:

  cfb:        Central fanbeam of the observation (int)
    OR
  utc_start:  UTC_START of the observation (str)
  


  Returns:
  --------
  fb:     The fan-beam number of the candidate
  '''
  
  obs_header = None
  if obs_header_file:
    obs_header = read_obs_header(obs_header_file = obs_header_file)

  elif utc_start:
    obs_header = read_obs_header(utc_start = utc_start)


  if utc_cand:
    utc_cand = DT.strptime(utc_cand, fmt_str(utc_cand))

  else:
    if utc_start is None or tcand is None:
      raise RuntimeError("Need either 'UTC_START AND tcand', or 'utc_cand'. ")

    utc_start = DT.strptime(utc_start, fmt_str(utc_start))
    if isinstance(tcand, float):
      tcand = TD(seconds = tcand)
    else:
      tcand = np.array([TD(seconds = i) for i in tcand])

    utc_cand = utc_start + tcand
  
  if fb_s is None:
    if obs_header is None:
      raise RuntimeError("Need to provide at least one of UTC_START or obs_header_file or fb_s")
    fb_s =  obs_header.FB_BEAM_SPACING 

  if cfb is None:
    if obs_header is None:
      raise RuntimeError("Need to provide at least one of UTC_START or obs_header_file or cfb")
    cfb = int((obs_header.NBEAM + 1)/2) + 1

  if isinstance(RA, str):
    RA = Angle(RA, unit=u.hourangle)
  if isinstance(DEC, str):
    DEC = Angle(DEC, unit=u.degree)
  
  
  if obs_header.DELAY_TRACKING:
    '''
    This peice of code will only work if the transiting center and the source are very nearby on sky
    I do not take into account the precession of sources, which should be roughly the same for nearby sources
    '''
    cfb_tracking_RA = Angle(obs_header.RA, unit=u.hourangle)
    cfb_tracking_DEC = Angle(obs_header.DEC, unit=u.degree)
    
    
    #delta_RA = RA * np.cos(DEC.radian) - cfb_tracking_RA * np.cos(cfb_tracking_DEC.radian)
    #delta_degrees = delta_RA.degree 
    #delta_fb = delta_degrees / fb_s
  
    MD_tracking_center = get_MD_from_RA_DEC(cfb_tracking_RA, cfb_tracking_DEC, utc_cand)
    MD_target = get_MD_from_RA_DEC(RA, DEC, utc_cand)

    delta_degrees = -1 *(MD_tracking_center - MD_target).deg
    delta_fb = delta_degrees / fb_s

  else:
    MD = get_MD_from_RA_DEC(RA, DEC, utc_cand)
    delta_fb = MD.deg / fb_s
  

  #Now since HA and MD are positive towards east, while FB numbers decrease towards east, we subtract delta_FB from the central fanbeam
  fb = cfb - delta_fb
  
  return fb





