#!/home/observer/miniconda2/bin/python

import ephem as e
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, FK5, ICRS
import astropy.units as u
from datetime import datetime as DT
from tqdm import trange, tqdm

#@profile
def work():
  fac = 0.2
  N_ra = int(360 * 5 * fac)
  N_dec = int(180 * 5 * fac)
  print "({},{})".format(N_ra, N_dec)
  
  ras = np.linspace(0., 360, N_ra)
  decs = np.linspace(-90., 90., N_dec)
  
  utc = DT.utcnow()
  
  shape = (N_ra, N_dec)

  era = np.empty(shape)
  edec = np.empty(shape)

  DECS_grid, RAS_grid = np.meshgrid(decs, ras)
  coords_2000 = SkyCoord(RAS_grid, DECS_grid, unit=(u.deg, u.deg), frame=FK5(equinox="J2000"))
  coords_now_aspy = coords_2000.transform_to(FK5(equinox=utc))
  #coords_now_aspy = coords_now_aspy.transform_to(ICRS)
  
  for i in trange(N_ra):
    for j in trange(N_dec):
      x = e.Equatorial(np.deg2rad(ras[i]), np.deg2rad(decs[j]), epoch=e.J2000)
      Ecoords_now = e.Equatorial(x, epoch=utc)
      
      era[i][j] = Ecoords_now.ra
      edec[i][j] = Ecoords_now.dec
      
  
  ra_diffs = -coords_now_aspy.ra.value + np.rad2deg(era)
  dec_diffs = -coords_now_aspy.dec.value + np.rad2deg(edec)
  tot_diffs = (ra_diffs**2 + dec_diffs**2)**0.5

  ra_diffs *= 3600
  dec_diffs *= 3600
  tot_diffs *= 3600

  cmap = 'seismic'
  extent = [decs[0], decs[-1],ras[0], ras[-1]]
  plt.figure()
  plt.imshow(ra_diffs, origin='lower', aspect='auto', cmap=cmap, extent=extent)
  plt.colorbar()
  plt.title("RA diffs")
  plt.xlabel("Declination")
  plt.ylabel("RA")
  
  
  plt.figure()
  plt.imshow(dec_diffs, origin='lower', aspect='auto', cmap='bwr', extent=extent)
  plt.colorbar()
  plt.title("DEC diffs")
  plt.xlabel("Declination")
  plt.ylabel("RA")
  
  
  plt.figure()
  plt.imshow(tot_diffs, origin='lower', aspect='auto', cmap=cmap, extent=extent)
  plt.colorbar()
  plt.title("Tot diffs")
  plt.xlabel("Declination")
  plt.ylabel("RA")
  plt.show()


work()
