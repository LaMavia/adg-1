# Using an edit-distance-like dynamic programming formulation, we can
# look for approximate occurrences of p in t.

import sys
import numpy as np
# from collections import deque
from array import array
from numpy._typing import NDArray
from tqdm import tqdm
from Bio import SeqIO
from sys import argv
import math

from contextlib import contextmanager
from time import perf_counter

@contextmanager
def catchtime(label: str):
    t1 = t2 = perf_counter()
    yield
    t2 = perf_counter()
    print(f'{label} took {t2 - t1:.3f} s')


def radixpass(a, b, r, n, k) :
  c = array("I", [0]*(k+1))
  for i in range(n) :
    c[r[a[i]]]+=1

  somme = 0
  for i in range(k+1):
    freq, c[i] = c[i], somme
    somme += freq

  for i in range(n) :
    b[c[r[a[i]]]] = a[i]
    c[r[a[i]]] += 1

def direct_kark_sort(s) :
  alphabet=[None, '$', 'A', 'C', 'G', 'T']
  k = len(alphabet)
  n = len(s)
  t = {None:0, '$': 1, 'A': 2, 'C': 3, 'G': 4, 'T': 5}
  SA = array('I', [0]*(n+3))
  kark_sort(array('I', [t[c] for c in s]+[0]*3), SA, n, k)
  return SA[:n]

def kark_sort(s, SA, n, K):
  n0  = (n+2) // 3
  n1  = (n+1) // 3
  n2  = n // 3
  n02 = n0 + n2
      
  SA12 = array('I', [0]*(n02+3))
  SA0  = array('I', [0]*n0)

  s12 = [i for i in range(n+(n0-n1)) if i%3 != 0] + [0,0,0] 
  s12 = array('I', s12)

  radixpass(s12, SA12, s[2:], n02, K)
  radixpass(SA12, s12, s[1:], n02, K)
  radixpass(s12, SA12, s, n02, K)

  name = 0
  c0, c1, c2 = -1, -1, -1
  for i in range(n02) :
    if s[SA12[i]] != c0 or s[SA12[i]+1] != c1 or s[SA12[i]+2] != c2 :
      name += 1
      c0 = s[SA12[i]]
      c1 = s[SA12[i]+1]
      c2 = s[SA12[i]+2]
    if SA12[i] % 3 == 1 :
      s12[SA12[i]//3] = name
    else :
      s12[SA12[i]//3 + n0] = name

  if name < n02 :
    kark_sort(s12, SA12, n02, name+1)
    for i in range(n02) :
      s12[SA12[i]] = i+1
  else :
    for i in range(n02) :
      SA12[s12[i]-1] = i

  s0 = array('I',[SA12[i]*3 for i in range(n02) if SA12[i]<n0])
  radixpass(s0, SA0, s, n0, K)
  
  p = j = k = 0
  t = n0 - n1
  while k < n :
    i = SA12[t]*3+1 if SA12[t]<n0 else (SA12[t] - n0)*3 + 2
    j = SA0[p] if p < n0 else 0

    if SA12[t] < n0 :
      test = (s12[SA12[t]+n0] <= s12[j//3]) if(s[i]==s[j]) else (s[i] < s[j])
    elif(s[i]==s[j]) :
      test = s12[SA12[t]-n0+1] <= s12[j//3 + n0] if(s[i+1]==s[j+1]) else s[i+1] < s[j+1]
    else :
      test = s[i] < s[j]

    if(test) :
      SA[k] = i
      t += 1
      if t == n02 :
        k += 1
        while p < n0 :
          SA[k] = SA0[p]
          p += 1
          k += 1
        
    else : 
      SA[k] = j
      p += 1
      if p == n0 :
        k += 1
        while t < n02 :
          SA[k] = (SA12[t] * 3) + 1 if SA12[t] < n0 else ((SA12[t] - n0) * 3) + 2
          t += 1
          k += 1
    k += 1

# Assume x is the string labeling rows of the matrix and y is the
# string labeling the columns

def suffixArray(t):
    with catchtime("suffixArray"):
        return direct_kark_sort(t)

def trace(D, x, y, yl: int, yr: int):
    ''' Backtrace edit-distance matrix D for strings x and y '''
    y_len = yr - yl
    i, j = len(x), y_len
    xscript = 0
    while i > 0:
        diag, vert, horz = sys.maxsize, sys.maxsize, sys.maxsize
        if i > 0 and j > 0:
            delt = 0 if x[i-1] == y[yl + j-1] else 1
            diag = D[i-1, j-1] + delt
        if i > 0:
            vert = D[i-1, j] + 1
        if j > 0:
            horz = D[i, j-1] + 1
        xscript += 1
        if diag <= vert and diag <= horz:
            # diagonal was best
            i -= 1; j -= 1
        elif vert <= horz:
            # vertical was best; this is an insertion in x w/r/t y
            i -= 1
        else:
            # horizontal was best
            j -= 1
    # j = offset of the first (leftmost) character of t involved in the
    # alignment
    return j, xscript

def allocDpArray(p: str, t_len: int):
    D = np.zeros((len(p)+1, t_len+1), dtype=np.uint32)
    # D = np.zeros((len(p)+1, t_len+1), dtype=int)
    D[1:, 0] = range(1, len(p)+1)

    return D

def kEditDp(p: str, t: str, tl: int, tr: int, k: int, D: NDArray):
    ''' Find and return the alignment of p to a substring of t with the
        fewest edits.  We return the edit distance, the offset of the
        substring aligned to, and the edit transcript.  If multiple
        alignments tie for best, we report the leftmost. '''
    t_len = tr - tl
    for i in range(1, len(p)+1):
        row_min = sys.maxsize
        for j in range(1, t_len+1):
            delt = 1 if p[i-1] != t[tl + j-1] else 0
            D[i, j] = (v := min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1))
            row_min = min(row_min, v)

        if row_min > k:
            return sys.maxsize, 0, 0

    # Find minimum edit distance in last row
    last_row = D[len(p), :]
    mnJ = last_row.argmin()
    mn = last_row[mnJ]
    # Backtrace; note: stops as soon as it gets to first row
    # t[lt:tr][:mnJ]
    off, xcript = trace(D, p, t, tl, tl + mnJ) #type: ignore
    # Return edit distance, offset into T, edit transcript
    return mn, off, xcript

def bwtFromSa(t, sa=None):
    ''' Given T, returns BWT(T) by way of the suffix array. '''
    bw = []
    dollarRow = None
    if sa is None:
        sa = suffixArray(t)
    for si in sa:
        if si == 0:
            dollarRow = len(bw)
            bw.append('$')
        else:
            bw.append(t[si-1])
    return ''.join(bw), dollarRow # return string-ized version of list bw

class FmCheckpoints(object):
    ''' Manages rank checkpoints and handles rank queries, which are
        O(1) time, with the checkpoints taking O(m) space, where m is
        length of text. '''
    
    def __init__(self, bw, cpIval=4):
        ''' Scan BWT, creating periodic checkpoints as we go '''
        self.cps = {}        # checkpoints
        self.cpIval = cpIval # spacing between checkpoints
        tally = {}           # tally so far
        # Create an entry in tally dictionary and checkpoint map for
        # each distinct character in text
        for c in bw:
            if c not in tally:
                tally[c] = 0
                self.cps[c] = []
        # Now build the checkpoints
        for i, c in enumerate(bw):
            tally[c] += 1 # up to *and including*
            if i % cpIval == 0:
                for c in tally.keys():
                    self.cps[c].append(tally[c])
    
    def rank(self, bw, c, row):
        ''' Return # c's there are in bw up to and including row '''
        if row < 0 or c not in self.cps:
            return 0
        i, nocc = row, 0
        # Always walk to left (up) when calculating rank
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1
        return self.cps[c][i // self.cpIval] + nocc

class FmIndex():
    ''' O(m) size FM Index, where checkpoints and suffix array samples are
        spaced O(1) elements apart.  Queries like count() and range() are
        O(n) where n is the length of the query.  Finding all k
        occurrences of a length-n query string takes O(n + k) time.
        
        Note: The spacings in the suffix array sample and checkpoints can
        be chosen differently to achieve different bounds. '''
    
    @staticmethod
    def downsampleSuffixArray(sa, n=4):
        ''' Take only the suffix-array entries for every nth suffix.  Keep
            suffixes at offsets 0, n, 2n, etc with respect to the text.
            Return map from the rows to their suffix-array values. '''
        ssa = {}
        for i, suf in enumerate(sa):
            # We could use i % n instead of sa[i] % n, but we lose the
            # constant-time guarantee for resolutions
            if suf % n == 0:
                ssa[i] = suf
        return ssa
    
    def __init__(self, t, cpIval=4, ssaIval=4):
        if t[-1] != '$':
            t += '$' # add dollar if not there already
        # Get BWT string and offset of $ within it
        sa = suffixArray(t)
        self.bwt, self.dollarRow = bwtFromSa(t, sa)
        # Get downsampled suffix array, taking every 1 out of 'ssaIval'
        # elements w/r/t T
        self.ssa = self.downsampleSuffixArray(sa, ssaIval)
        self.slen = len(self.bwt)
        # Make rank checkpoints
        self.cps = FmCheckpoints(self.bwt, cpIval)
        # Calculate # occurrences of each character
        tots = dict()
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1
        # Calculate concise representation of first column
        self.first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            self.first[c] = totc
            totc += count
    
    def count(self, c):
        ''' Return number of occurrences of characters < c '''
        if c not in self.first:
            # (Unusual) case where c does not occur in text
            for cc in sorted(self.first.iterkeys()):
                if c < cc: return self.first[cc]
            return self.first[cc]
        else:
            return self.first[c]
    
    def range(self, p):
        ''' Return range of BWM rows having p as a prefix '''
        l, r = 0, self.slen - 1 # closed (inclusive) interval
        for i in range(len(p)-1, -1, -1): # from right to left
            l = self.cps.rank(self.bwt, p[i], l-1) + self.count(p[i])
            r = self.cps.rank(self.bwt, p[i], r)   + self.count(p[i]) - 1
            if r < l:
                break
        return l, r+1
    
    def resolve(self, row):
        ''' Given BWM row, return its offset w/r/t T '''
        # with catchtime(f'resolve {row}'):
        def stepLeft(row):
            ''' Step left according to character in given BWT row '''
            c = self.bwt[row]
            return self.cps.rank(self.bwt, c, row-1) + self.count(c)
        nsteps = 0
        while row not in self.ssa:
            row = stepLeft(row)
            nsteps += 1
        return self.ssa[row] + nsteps
    
    def hasSubstring(self, p):
        ''' Return true if and only if p is substring of indexed text '''
        l, r = self.range(p)
        return r > l
    
    def hasSuffix(self, p):
        ''' Return true if and only if p is suffix of indexed text '''
        l, r = self.range(p)
        off = self.resolve(l)
        return r > l and off + len(p) == self.slen-1
    
    def occurrences(self, p):
        ''' Return offsets for all occurrences of p, in no particular order '''
        l, r = self.range(p)
        # print(f"[FMI.occurrences] range_len={r - l + 1}")
        return (self.resolve(x) for x in range(l, r))

# First we make a function that splits a string p up into a set of
# non-overlapping, non-empty substrings.

def partition(p, pieces=2):
    assert len(p) >= pieces
    base, mod = len(p) // pieces, len(p) % pieces
    idx = 0
    ps = []
    modAdjust = 1
    for i in range(0, pieces):
        if i >= mod: modAdjust = 0
        newIdx = idx + base + modAdjust
        ps.append((p[idx:newIdx], idx))
        idx = newIdx
    return ps

B = 0
NB = 0

def queryIndexEdit(p, t, k, index):
    ''' Look for occurrences of p in t with up to k edits using an
        index combined with dynamic-programming alignment. '''
    global B, NB
    def aux(D: NDArray | None, dp_range: int, lf: int, rt: int):
        # left edge of T to include in DP matrix

        if (hit_range := rt - lf) > dp_range or D is None:
            D = allocDpArray(p, hit_range)
            dp_range = hit_range
        # with catchtime(f"[kEditDp] range={rt-lf + 1}"):
        mn, off, xcript = kEditDp(p, t, lf,rt, k, D)
        off += lf
        ret = None
        if mn <= k:
            ret = (mn, off, off + xcript)
        return D, lf, rt, dp_range, ret


    for part, poff in partition(p, k+1):
        D : NDArray | None = None
        dp_range = 0
        last_lf = last_rt = lf = rt = 0
        for hit in index.occurrences(part): # query index w/ partition
            lf = max(0, hit - poff - k)
            # right edge of T to include in DP matrix
            rt = min(len(t), hit - poff + len(p) + k)
            if lf < last_rt:
                B += 1
                continue

            NB += 1
            last_rt = rt
            last_lf = lf
            D, lf, rt, dp_range, ret = aux(D, dp_range, lf, rt)
            if ret is not None:
                yield ret

        if last_rt != rt or last_lf != lf:
            D, lf, rt, dp_range, ret = aux(D, dp_range, lf, rt)
            if ret is not None:
                yield ret
                
def main():
    global B, NB
    seq_rec=next(SeqIO.parse(argv[1], "fasta"))
    t = str(seq_rec.seq)
    with catchtime("fm-index"):
        index = FmIndex(t, 11, 11)

    fout = open(argv[3], "w")
    reads = list(SeqIO.parse(argv[2], "fasta"))
    bar = tqdm(reads)
    for read in bar:
        for hit in queryIndexEdit(str(read.seq), t, len(read.seq)//9, index):
            fout.write("{}\t{}\t{}\n".format(read.id, hit[1], hit[2]))
            break
        bar.set_postfix({'bail': f'{B/(NB + B) * 100:.2f}%'})
    fout.close()

if __name__ == "__main__":
    main()
