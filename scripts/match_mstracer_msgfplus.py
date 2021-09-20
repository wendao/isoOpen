import sys

lines = open(sys.argv[2], 'r').readlines()
t = {}
for n, tag in enumerate(lines[0].strip().split('\t')):
    t[tag] = n
#{'SpecID': 1, 'DeNovoScore': 10, 'EValue': 13, 'Precursor': 4, 'ScanNum': 2, 'SpecEValue': 12, 'IsotopeError': 5, 'PrecursorError(ppm)': 6, 'FragMethod': 3, 'PepQValue': 15, 'Charge': 7, 'MSGFScore': 11, 'QValue': 14, 'Peptide': 8, 'Protein': 9, '#SpecFile': 0}
for l in lines[1:]:
    es = l.strip().split('\t')
    ev = float(es[t["EValue"]])
    if ev < 0.01:




lines = open(sys.argv[1], 'r').readlines()
t = {}
for n, tag in enumerate(lines[0].strip().split('\t')):
    t[tag] = n
#{'rtApex': 4, 'massCalib': 11, 'rtStart': 3, 'averagineCorr': 9, 'mostAbundantMz': 1, 'charge': 2, 'rtEnd': 5, 'nIsotopes': 7, 'intensitySum': 13, 'mass': 10, 'nScans': 8, 'fwhm': 6, 'intensityApex': 12, 'mz': 0}
for l in lines[1:]:
    es = l.strip().split('\t')
    rtStart = float(es[t["rtStart"]])
    rtEnd = float(es[t["rtEnd"]])
    rtApex = float(es[t["rtApex"]])
    print rtStart, rtEnd
