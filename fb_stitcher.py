
import numpy as np
from sigpyproc.Readers import FilReader as F
import argparse, os, sys, glob
from tqdm import tqdm

from get_fanbeam_from_coord import get_fanbeam_from_coordinates as G

def main():

  print("Getting the observation header information")
  utc_start = args.utc_start
  last_beam = args.nbeams
  any_fil = glob.glob(os.path.join(args.basedir, "BEAM_*/{}.fil".format(utc_start)) )
  if len(any_fil) == 0:
    raise RuntimeError("Could not find any filterbanks in the basedir")

  f = F(any_fil[0])

  total_tobs = f.header.nsamples * f.header.tsamp
  res = int(args.res / f.header.tsamp) * f.header.tsamp   #Making the tstep resolution an integral multiple of the sampling time, so that we do not need to step by fractional samples

  tt = np.arange(0, total_tobs, res)

  source_name = args.sname
  source_ra = args.ra.strip()
  source_dec = args.dec.strip()

  if args.load_track:
    print("Reading which beams the source would have been in during the observation from {}".format(args.load_track))
    source_track = np.loadtxt(args.load_track)
    source_beams = source_track[:,0]
    source_times = source_track[:,1]

    print("Using the time resolution from the track file instead of the value specified as argument (if any)")
    res = int( ( source_times[1] - source_times[0] ) / f.header.tsamp ) * f.header.tsamp
  else:
    print("Getting which beams the source would have been in during the observation")
    source_beams =  G(source_ra, source_dec, utc_start=utc_start, tcand=tt)
    mask = (2 <= source_beams) & (source_beams <= last_beam)

    source_beams = source_beams[mask]
    source_times = tt[mask]

  if len(source_beams) ==  0:
    print("The source did not transit any fan beam in this observation")
    exit()

  start_samp = int(source_times[0] / f.header.tsamp)
  nsamp = int(res / f.header.tsamp)

  update_hdr = {
      'tstart': f.header.mjdAfterNsamps(start_samp),
      'source_name': source_name,
      'ra':         source_ra,
      'dec':        source_dec
      }
  of = f.header.prepOutfile("{0}_{1}_{2}_{3}.fil".format(source_name, source_ra, source_dec, utc_start), updates=update_hdr)


  for ii, tstep in enumerate(tqdm(source_times)):
    fbeam = source_beams[ii]
    lbeam = int(fbeam)
    rbeam = lbeam + 1

    ss = int(tstep / f.header.tsamp)

    beam_frac = fbeam - lbeam

    lbf = os.path.join(args.basedir + "BEAM_{0:03g}/{1}.fil".format(lbeam, utc_start))
    rbf = os.path.join(args.basedir + "BEAM_{0:03g}/{1}.fil".format(rbeam, utc_start))

    lf, rf = None, None
    ld, rd = 0, 0
    if os.path.exists(lbf):
      lf = F(lbf)
      ld = lf.readBlock(ss, nsamp)
    else:
      if lbeam > 1:
        print("Could not find BEAM_{0:03g} -- skipping".format(lbeam))

    if os.path.exists(rbf):
      rf = F(rbf)
      rd = rf.readBlock(ss,  nsamp)
    else:
      if rbeam < last_beam:
        print("Could not find BEAM_{0:03g} -- skipping".format(rbeam))

    if lf and rf:
      od = ld * (1-beam_frac) + rd * beam_frac
      od_mean = np.mean(od, axis=1)
      od -= od_mean[:, None]
      norm = np.sqrt( (1-beam_frac)**2 + beam_frac**2  )
      od /= norm
      od += od_mean[:, None]

    elif lf:
      od = ld
    elif rf:
      od = rd
    else:
      print("Something went terribly wrong")

    of.cwrite(od.T.ravel().astype(f.header.dtype, casting='unsafe'))

  of.close()
  if args.save_track:
    outname = "Track_{0}_{1}_{2}_{3}.txt".format(source_name, source_ra, source_dec, utc_start)
    print("Saving the source track info to {}".format(outname))
    np.savetxt(outname, np.hstack([source_beams[:,None], source_times[:,None]]), fmt='%.3f', header="Beams\tTimes")
  print("Done")

if __name__ == '__main__':
  for i, arg in enumerate(sys.argv):
    if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
  a = argparse.ArgumentParser()
  a.add_argument("-utc_start", type=str, help="UTC_START of the observation", required=True)
  a.add_argument("-basedir", type=str, help="Base directory which contains all the beams (like /data/mopsr/archives/<utc>/FB/)", required=True)
  a.add_argument("-sname", type=str, help="Source name to stitch for", default="Stitching_test_source")
  a.add_argument("-ra", type=str, help="RA of the source in HH:MM:SS.ff format", required=True)
  a.add_argument("-dec", type=str, help="DEC of the source in DD:MM:SS.ff format", required=True)
  a.add_argument("-res", type=float, help="Time resolution (in seconds) at which to sample the fan-beams (def=0.1)", default=0.1)
  a.add_argument("-nbeams", type=int, help="Number of fanbeams in the observation (def=352)", default=352)
  a.add_argument("-load_track", type=str, help="Load the source track from a file instead of calculating -- takes the path to the track file as argument", default=None)
  a.add_argument("-save_track", action='store_true', help="Invoke this option if you would like to save the source track in a file on disk", default=False)
  args = a.parse_args()
  main()


