# Using an edit-distance-like dynamic programming formulation, we can
# look for approximate occurrences of p in t.

import sys
import numpy
from collections import deque

# Assume x is the string labeling rows of the matrix and y is the
# string labeling the columns

def trace(D, x, y, k):
    ''' Backtrace edit-distance matrix D for strings x and y '''
    # idziemy od komórek ostatniego wiersza, które są <=k, po czym robimy backtracking,
    # ale musimy trzymać cały zbiór rozwiązań, ponieważ w niektórych przypadkach możemy
    # mieć alternatywne ścieżki

    done_positions = []
    positions = deque([(len(x), j, []) for j in range(len(y)) if j <= k])

    while len(positions) > 0:
        i, j, xscript = positions.pop()

        while i > 0:
            diag, vert, horz = sys.maxsize, sys.maxsize, sys.maxsize
            delt = None
            if i > 0 and j > 0:
                delt = 0 if x[i-1] == y[j-1] else 1
                diag = D[i-1, j-1] + delt
            if i > 0:
                vert = D[i-1, j] + 1
            if j > 0:
                horz = D[i, j-1] + 1
            
            min_val = min(diag, vert, horz)
            if diag == min_val:
                positions.append((i - 1, j - 1, [*xscript, 'R' if delt == 1 else 'M']))
            if vert == min_val:
                positions.append((i-1, j, [*xscript, 'I']))
            if horz == min_val:
                positions.append((i, j-1, [*xscript, 'D']))
    # j = offset of the first (leftmost) character of t involved in the
    # alignment
    return [(j, (''.join(xscript))[::-1]) for _, j, xscript in positions]

def kEditDp(p, t, k: int):
    ''' Find and return the alignment of p to a substring of t with the
        fewest edits.  We return the edit distance, the offset of the
        substring aligned to, and the edit transcript.  If multiple
        alignments tie for best, we report the leftmost. '''
    D = numpy.zeros((len(p)+1, len(t)+1), dtype=int)
    # Note: First row gets zeros.  First column initialized as usual.
    D[1:, 0] = range(1, len(p)+1)
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            delt = 1 if p[i-1] != t[j-1] else 0
            D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
    # Find minimum edit distance in last row
    mnJ, mn = None, len(p) + len(t)
    for j in range(len(t)+1):
        if D[len(p), j] < mn:
            mnJ, mn = j, D[len(p), j]
    # Backtrace; note: stops as soon as it gets to first row
    results = trace(D, p, t[:mnJ], k)
    # Return edit distance, offset into T, edit transcript
    return mn, results, D

def suffixArray(s):
    ''' Given T return suffix array SA(T).  Uses "sorted"
        function for simplicity, which is probably very slow. '''
    satups = sorted([(s[i:], i) for i in range(len(s))])
    return list(map(lambda x: x[1], satups)) # extract, return just offsets

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
            for cc in sorted(self.first.keys()):
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
        return [ self.resolve(x) for x in range(l, r) ]

def partition(p, k):
  return [p[i:i+k] for i in range(len(p), k)]

# def queryIndexEdit(p, t, k, index):
#     ''' Look for occurrences of p in t with up to k edits using an
#         index combined with dynamic-programming alignment. '''
#     l = index.ln
#     occurrences = []
#     seen = set()     # for avoiding reporting same hit twice
#     for part, poff in partition(p, k+1):
#         for hit in index.occurrences(part): # query index w/ partition
#             # left edge of T to include in DP matrix
#             lf = max(0, hit - poff - k)
#             # right edge of T to include in DP matrix
#             rt = min(len(t), hit - poff + len(p) + k)
#             mn, off, xcript = kEditDp(p, t[lf:rt])
#             off += lf
#             if mn <= k and (mn, off) not in seen:
#                 occurrences.append((mn, off, xcript))
#                 seen.add((mn, off))
#     return occurrences

