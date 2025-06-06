import sys
import argparse
import numpy as np
from math import fabs
from collections import defaultdict

hits_msf = []
mark_msf = {}
#fit Rt and scan, Rt in second
lines = open("psm.tsv", 'r').readlines()
t = {}
for n, tag in enumerate(lines[0].strip().split('\t')):
    t[tag] = n
#{'Next AA': 5, 'Protein Description': 32, 'Protein Start': 22, 'Charge': 7, 'Spectrum File': 1, 'Delta Mass': 15, 'Number of Enzymatic Termini': 20, 'Peptide Length': 6, 'Assigned Modifications': 25, 'Peptide': 2, 'Observed M/Z': 11, 'Protein ID': 29, 'Calculated M/Z': 14, 'Expectation': 16, 'Hyperscore': 17, 'PeptideProphet Probability': 19, 'Calculated Peptide Mass': 13, 'Mapped Proteins': 34, 'Observed Mass': 9, 'Calibrated Observed M/Z': 12, 'Observed Modifications': 26, 'Calibrated Observed Mass': 10, 'Modified Peptide': 3, 'Mapped Genes': 33, 'Entry Name': 30, 'Gene': 31, 'Is Unique': 27, 'Prev AA': 4, 'Protein End': 23, 'Nextscore': 18, 'Spectrum': 0, 'Intensity': 24, 'Number of Missed Cleavages': 21, 'Protein': 28, 'Retention': 8}
rs = []
for l in lines[1:]:
    es = l.strip().split('\t')
    mz = float(es[t["Observed M/Z"]])
    rt = float(es[t["Retention"]])
    spectrum = es[t["Spectrum"]]
    scan = int(spectrum.split('.')[1])
    charge = int(es[t["Charge"]])
    hits_msf.append((scan, mz, rt, charge))
    if scan % 10 == 0:
        rs.append(rt/scan)
r = np.median(rs)

hits_msg = []
mark_msg = {}
lines = open(sys.argv[2], 'r').readlines()
t = {}
for n, tag in enumerate(lines[0].strip().split('\t')):
    t[tag] = n
#{'SpecID': 1, 'DeNovoScore': 10, 'EValue': 13, 'Precursor': 4, 'ScanNum': 2, 'SpecEValue': 12, 'IsotopeError': 5, 'PrecursorError(ppm)': 6, 'FragMethod': 3, 'PepQValue': 15, 'Charge': 7, 'MSGFScore': 11, 'QValue': 14, 'Peptide': 8, 'Protein': 9, '#SpecFile': 0}
for l in lines[1:]:
    es = l.strip().split('\t')
    ev = float(es[t["EValue"]])
    if ev < 0.01:
        mz = float(es[t["Precursor"]])
        charge = int(es[t["Charge"]])
        scan = es[t["ScanNum"]]
        rt = r*int(scan)
        hits_msg.append((scan, mz, rt, charge))

#match, rt in min
ppm = 1e-6
rt_cut = 30.0

parser = argparse.ArgumentParser(description='Match MSTracer features with MSGF+ results')
parser.add_argument('mstracer_file', help='MSTracer feature file')
parser.add_argument('msgfplus_file', help='MSGF+ results file')
args = parser.parse_args()

print("Matching")
lines = open(args.mstracer_file, 'r').readlines()
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
    for scan_, mz_, rt_, c_ in hits_msf:
        if charge != c_: continue
        if fabs(mz_-mz) > 10.0*ppm*mz_:
            if fabs(mz_-mzI) > 10.0*ppm*mz_: continue
        #print rt_, rtStart, rtEnd
        if rt_<rtStart-rt_cut or rt_>rtEnd+rt_cut: continue
        print("Hit MSF:", scan_, mz_, rtStart, rtEnd, rt_)
        mark_msf[scan_] = 1
    #from msg
    for scan_, mz_, rt_, c_ in hits_msg:
        if charge != c_: continue
        if fabs(mz_-mz) > 10.0*ppm*mz_:
            if fabs(mz_-mzI) > 10.0*ppm*mz_: continue
        #print rt_, rtStart, rtEnd
        if rt_<rtStart-rt_cut or rt_>rtEnd+rt_cut: continue
        print("Hit MSG:", scan_, mz_, rtStart, rtEnd, rt_)
        mark_msg[scan_] = 1

t1 = len(mark_msf.keys())
t2 = len(hits_msf)
print(t1, t2, float(t1)/t2)
t1 = len(mark_msg.keys())
t2 = len(hits_msg)
print(t1, t2, float(t1)/t2)
