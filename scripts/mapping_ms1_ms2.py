from pyteomics import mzml
from collections import defaultdict
import os, sys

mzml_file_path = sys.argv[1]

ms1_to_ms2 = []

number_ms1 = 0
with mzml.read(mzml_file_path) as reader:
    for idx, spectrum in enumerate(reader):
        ms_level = spectrum.get('ms level')
        scan_number = spectrum.get('id', f"index={idx}").split('=')[-1]
        if ms_level == 1:
            number_ms1 += 1
            if number_ms1>1: print()
            print(number_ms1, end="\t")
            print(scan_number, end="\t")
            flag = ""
        elif ms_level == 2:
            print(flag+str(scan_number), end="")
            flag = ","
print()
