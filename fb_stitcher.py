
import numpy as np
from sigpyproc.Readers import FilReader as F
import argparse, os

def main():
  f = F(beam1)
  total_tobs = f.header.length

  tt = np.arange(0, total_tobs, args.res)

