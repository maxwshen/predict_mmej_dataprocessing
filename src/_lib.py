# Stores project-specific functions

import _config
import re
import numpy as np

test = 1

def translate(loc):
  # Input: chr15:40987440-40987463
  # Output: chro, start, end
  w = loc.split(':')
  chro = w[0].replace('chr', '')
  start = w[1].split('-')[0]
  end = w[1].split('-')[1]
  return chro, int(start), int(end)


def exp_fold_name(nm):
  # Input: SRR3702721
  # Output: Folder (SRR3702)
  if len(nm) == 10:
    return nm[:7]
  else:
    print '\tERROR:', nm, 'not valid (not 10 chars long)'
    return ''

def find_matching_spacer(loc, spacers):
  #
  chro, start, end = translate(loc)
  fnd = []
  for i in range(len(spacers)):
    if spacers[i]['chr'] == chro:
      if spacers[i]['start'] == start:
        fnd = [i] + fnd
      if spacers[i]['start'] - 3 <= start <= spacers[i]['start'] + 3:
        fnd.append(i)
  if len(fnd) == 0:
    print '\tNo matching spacer found to', loc
    return None
  else:
    return spacers[fnd[0]]

def read_cigarfile(fn):
  data = {}
  with open(fn) as f:
    for i, line in enumerate(f):
      if line.startswith('>'):
        curr_count = line.strip('>')
      else:
        curr_cigar = line.strip()
        data[curr_cigar] = curr_count
  return data

def split_cigar(cigar):
  types = []
  for i in range(len(cigar)):
    if cigar[i].isalpha():
      types.append(cigar[i])
  lens = re.split(r'M|D|I|N|S|H|P|=|X', cigar)
  lens = [int(s) for s in lens[:-1]]  # cigars end in letters, resulting in final extra empty string
  # print types, lens
  return types, lens

def convert_dict_to_list(d):
  l = [0] * (max(d.keys()) + 1)
  for i in range(len(l)):
    if i in d:
      l[i] += d[i]
  return np.array(l)