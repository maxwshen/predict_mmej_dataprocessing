#

import _config, _lib
import sys, os, fnmatch, datetime, subprocess, pickle, re
import numpy as np
from collections import defaultdict
from mylib import util
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'a_get_cigar/' + _config.NM + '/'
NAME = util.get_fn(__file__)

# Functions

def convert_max(data):
  dels = {}
  for key in data:
    types, lens = _lib.split_cigar(key)
    if 'D' not in types:
      continue 

    # find deletion length
    del_lens = []
    for i in range(len(types)):
      if types[i] == 'D':
        del_lens.append(lens[i])
    del_len = max(del_lens)

    # summarize into dict
    if del_len not in dels:
      dels[del_len] = 0
    dels[del_len] += 1
  return _lib.convert_dict_to_list(dels)

def convert_first(data):
  dels = {}
  for key in data:
    types, lens = _lib.split_cigar(key)
    if 'D' not in types:
      continue 

    # find deletion length
    del_lens = []
    for i in range(len(types)):
      if types[i] == 'D':
        del_lens.append(lens[i])
    del_len = del_lens[0]

    # summarize into dict
    if del_len not in dels:
      dels[del_len] = 0
    dels[del_len] += 1
  return _lib.convert_dict_to_list(dels)

def convert_last(data):
  dels = {}
  for key in data:
    types, lens = _lib.split_cigar(key)
    if 'D' not in types:
      continue 

    # find deletion length
    del_lens = []
    for i in range(len(types)):
      if types[i] == 'D':
        del_lens.append(lens[i])
    del_len = del_lens[-1]

    # summarize into dict
    if del_len not in dels:
      dels[del_len] = 0
    dels[del_len] += 1
  return _lib.convert_dict_to_list(dels)

def convert_max_fix_artifacts(data):
  dels = {}
  num_artifacts_fixed = 0
  for key in data:
    types, lens = _lib.split_cigar(key)
    if 'D' not in types:
      continue 

    # find deletion length
    del_lens = []
    idx = []
    for i in range(len(types)):
      if types[i] == 'D':
        del_lens.append(lens[i])
        idx.append(i)
    del_len = max(del_lens)

    whichD = idx[del_lens.index(max(del_lens))]

    # check left side
    if whichD > 1:
      if types[whichD - 2] == 'D' and types[whichD - 1] == 'M':
        if lens[whichD - 1] <= 3:
          del_len += lens[whichD - 2]
          # print 'adding', lens[whichD - 2], 'for', types, lens
          num_artifacts_fixed += 1

    # check right side
    if len(types) - whichD > 2:
      if types[whichD + 2] == 'D' and types[whichD + 1] == 'M':
        if lens[whichD + 1] <= 3:
          del_len += lens[whichD + 2]
          # print 'adding', lens[whichD + 2], 'for', types, lens
          num_artifacts_fixed += 1

    # summarize into dict
    if del_len not in dels:
      dels[del_len] = 0
    dels[del_len] += 1
  print '\tFixed', num_artifacts_fixed, 'artifacts'
  return _lib.convert_dict_to_list(dels)

def convert_spacers(inp_dir, out_dir, spacers):
  ctr = 0
  aad = []
  for fn in os.listdir(inp_dir):
    if '.txt' in fn:
      data = _lib.read_cigarfile(inp_dir + fn)
      
      nms, all_dels = [], []

      dels = convert_max(data)
      nms.append('max')
      all_dels.append(dels)
      
      dels = convert_max_fix_artifacts(data)
      nms.append('max | fixed artifacts')
      all_dels.append(dels)

      dels = convert_first(data)
      nms.append('first')
      all_dels.append(dels)

      dels = convert_last(data)
      nms.append('last')
      all_dels.append(dels)

      ymax = max([max(s) for s in all_dels])
      xmax = max([len(s) for s in all_dels])

      with PdfPages(out_dir + 'profile_' + fn.strip('.txt').replace(':', '_') + '.pdf') as pdf:
        for s in range(len(all_dels)):
          dels = all_dels[s]
          plt.bar(range(0, len(dels)), dels, align = 'center', color = '#D00000')
          plt.xlim([-2, xmax])
          plt.ylim([0, ymax + 10])
          plt.ylabel(nms[s])
          pdf.savefig()
          plt.close()

      aad.append(all_dels)
      
      ctr += 1
      if ctr == 3:
        for j in range(1, 31):
          print 'max j', np.std([aad[0][0][j], aad[1][0][j], aad[2][0][j] ])
          print 'maxf j', np.std([aad[0][1][j], aad[1][1][j], aad[2][1][j] ])
          print 'first j', np.std([aad[0][2][j], aad[1][2][j], aad[2][2][j] ])
          print 'last j', np.std([aad[0][3][j], aad[1][3][j], aad[2][3][j] ])
      if ctr > 5:
        return

  return


@util.time_dec
def main(inp_dir, out_dir, run = True):
  print NAME  
  util.ensure_dir_exists(out_dir)
  if not run:
    print '\tskipped'
    return out_dir

  # Function calls
  with open(inp_dir + 'spacers.pkl') as f:
    spacers = pickle.load(f)

  convert_spacers(inp_dir, out_dir, spacers)

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
