import sys
import numpy as np
from math import fabs, log
from collections import defaultdict

hits = []
mark = {}
lines = open(sys.argv[2], 'r').readlines()
t = {}
for n, tag in enumerate(lines[0].strip().split('\t')):
    t[tag] = n
#{'R2.20181026_1TO1': 12, 'index': 0, 'description': 2, 'NP.20181026_1TO1': 11, 'sequence': 4, 'symbol': 3, 'ipi': 1, 'INT.20181026_1TO1': 10, 'charge': 6, 'mass': 5, 'entry': 13, 'link': 15, 'IR.20181026_1TO1': 8, 'segment': 7, 'LR.20181026_1TO1': 9, 'RT': 14}
for l in lines[1:]:
    es = l.strip().split('\t')
    mass = float(es[t["mass"]])
    charge = int(es[t["charge"]])
    rt = float(es[t["RT"]])*60.0
    mz = (mass+1.0)/charge
    r = float(es[t["IR.20181026_1TO1"]])
    scan = int(es[t["index"]])
    try:
        if fabs(log(r))<0.693: #0.5 ~ 2
            hits.append((scan,mz,rt,charge))
    except:
        pass

#match, rt in min
ppm = 1e-6
rt_cut = 30.0

print "Matching"
lines = open(sys.argv[1], 'r').readlines()
t = {}
for n, tag in enumerate(lines[0].strip().split('\t')):
    t[tag] = n
#{'rtApex': 4, 'massCalib': 11, 'rtStart': 3, 'averagineCorr': 9, 'mostAbundantMz': 1, 'charge': 2, 'rtEnd': 5, 'nIsotopes': 7, 'intensitySum': 13, 'mass': 10, 'nScans': 8, 'fwhm': 6, 'intensityApex': 12, 'mz': 0}
for l in lines[1:]:
    es = l.strip().split('\t')
    mz = float(es[t["mz"]])
    mzI = float(es[t["mostAbundantMz"]])
    rtStart = float(es[t["rtStart"]])*60.0
    rtEnd = float(es[t["rtEnd"]])*60.0
    rtApex = float(es[t["rtApex"]])*60.0
    charge = int(es[t["charge"]])
    #from msf
    for scan_, mz_, rt_, c_ in hits:
        if charge != c_: continue
        if fabs(mz_-mz) > 20.0*ppm*mz_:
            if fabs(mz_-mzI) > 20.0*ppm*mz_: continue
        #print rt_, rtStart, rtEnd
        if rt_<rtStart-rt_cut or rt_>rtEnd+rt_cut: continue
        print "Hit DTA:", scan_, mz_, rtStart, rtEnd, rt_
        mark[scan_] = 1

t1 = len(mark.keys())
t2 = len(hits)
print t1, t2, float(t1)/t2
