#

import _config, _lib
import sys, os, fnmatch, datetime, subprocess, pickle
import numpy as np, pandas as pd
from collections import defaultdict
from mylib import util


# Default params
DEFAULT_INP_DIR = _config.DATA_DIR
NAME = util.get_fn(__file__)

# Functions

def load_spacers(inp_dir):
  spacers = defaultdict(dict)
  df = pd.read_csv(inp_dir + _config.d.SPACER_FN)

  print '\tFound', df.shape[0], 'spacers'
  for i in range(df.shape[0]):
    num = df[_config.d.SPC_NUM_COL][i]
    loc = df[_config.d.SPC_LOC_COL][i]
    chro, start, end = _lib.translate(loc)
    spacers[i]['chr'] = chro
    spacers[i]['start'] = int(start)
    spacers[i]['end'] = int(end)
    spacers[i]['num'] = int(num)
  return spacers

def find_exps(inp_dir, spacers):
  df = pd.read_csv(inp_dir + _config.d.INFO_FN)

  print '\tUsing nm regex:', _config.NM
  
  for spc in spacers.values():
    spc['runs'] = []
    spc['libnms'] = []

  for i in range(df.shape[0]):
    lib_nm = df[_config.d.INF_LIB_COL][i]
    run_nm = df[_config.d.INF_RUN_COL][i]

    if _config.NM in lib_nm:
      loc = lib_nm.split('_')[1]
      spc = _lib.find_matching_spacer(loc, spacers)
      if spc is not None:
        spc['runs'].append(run_nm)
        spc['libnms'].append(lib_nm)

  for spc in spacers.values():
    print '\t\t', spc['num'], ':', len(spc['runs']),
    if len(spc['runs']) == 0:
      print ' *'
    else:
      print ''
  return spacers

def get_cigars(inp_dir, out_dir, spacers):
  print '\tGetting cigars...'
  timer = util.Timer(total = len(spacers))
  for spc in spacers.values():
    for i in range(len(spc['runs'])):
      cigars = {}
      run = spc['runs'][i]
      foldnm = _lib.exp_fold_name(run)
      fn = _config.DATA_DIR + foldnm + '/' + run + '.sam'
      num_aligns, num_kept = 0.0, 0.0
      with open(fn) as f:
        for _, line in enumerate(f):
          if not line.startswith('@'):
            num_aligns += 1
            chro, start = line.split()[2], int(line.split()[3])
            cigar = line.split()[5]

            if chro == 'chr' + spc['chr']:
              if spc['start']-300 <= start <= spc['start']+300:
                if cigar not in cigars:
                  cigars[cigar] = 0
                cigars[cigar] += 1
                num_kept += 1

      frac_kept = num_kept / num_aligns
      if frac_kept < 0.80 and num_aligns > 1000:
        print '\tWARNING: Kept:', num_kept / num_aligns, 'of', num_aligns, 'for spacer', spc['num'], ':', run, spc['libnms'][i]

      out_fn = out_dir + spc['libnms'][i] + '.txt'
      with open(out_fn, 'w') as f:
        for cigar in cigars:
          f.write('>' + str(cigars[cigar]) + '\n' + cigar + '\n')
    timer.update()

  return


@util.time_dec
def main(inp_dir, out_dir, run = True):
  print NAME  
  util.ensure_dir_exists(out_dir)
  if not run:
    print '\tskipped'
    return out_dir

  # Function calls
  spacers = load_spacers(inp_dir)
  spacers = find_exps(inp_dir, spacers)
  get_cigars(inp_dir, out_dir, spacers)
  
  print '\tWriting pickle...'
  with open(out_dir + 'spacers.pkl', 'wb') as f:
    pickle.dump(spacers, f)

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/' + _config.NM + '/')