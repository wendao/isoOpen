import sys
import argparse
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit
from spectutilio import *
import numpy as np
from math import fabs
from spectrum_utils import fragment_annotation as fa, proforma
from ms2pip.single_prediction import SinglePrediction

#hyper parameters
alpha = 2000.0  #sharpness
bin_size = 0.01 #Da
nbin = int((max_mz - min_mz)/bin_size)
halfL = int(nbin/2)
print( "bin_size=", bin_size, "nbin=", nbin )
xbin = np.linspace(min_mz, max_mz, nbin+1)
x_mid = (min_mz+max_mz)/2

def gaussian(x, alpha, r):
  return np.exp(-alpha*np.power((x - r), 2.))

def parabola(x, a, b, c):
    return a * (x - b)**2 + c

def remove_ions_from(ori_mz, pat_mz, eps=0.008):
    result = []
    for element in ori_mz:
        dist0 = [ fabs(element - pattern) for pattern in pat_mz ]
        dist1 = [ fabs(element - pattern - 1.003355) for pattern in pat_mz ]
        dist2 = [ fabs(element - pattern - 1.003355*2) for pattern in pat_mz ]
        if np.min(dist0)>eps and np.min(dist1)>eps and np.min(dist2)>eps: result.append(element)
    return result

def extract_intensity_from(ori_mz, ori_int, pat_mz, eps=0.005, fill=0.0):
    pat_int = []
    for pattern in pat_mz:
        matched_index = None
        for i, mz in enumerate(ori_mz):
            if abs(mz - pattern) < eps:
                matched_index = i
                break
        # If a matching element is found, append its intensity to pat_int; otherwise, append 0
        if matched_index is not None:
            pat_int.append(ori_int[matched_index])
        else:
            pat_int.append(fill)
    return pat_int

def fit_mz_by_sig_peak( sig, p_loc, N=3 ):
    global bin_size
    x = [ i*bin_size for i in range(-N, N+1) ]
    y = sig[ p_loc-N : p_loc+N+1 ]
    #print(x, y)
    initial_guess = [-1000, 0, 10]
    params, covariance = curve_fit(parabola, x, y, p0=initial_guess)
    a_fit, b_fit, c_fit = params
    #print( a_fit, b_fit, c_fit )
    acc_mz = b_fit

    #draw qc
    x_fit = np.linspace(min(x)-bin_size, max(x)+bin_size, 100*N)
    y_fit = parabola(x_fit, a_fit, b_fit, c_fit)

    plt.scatter(x, y, label='corr score', color='blue', s=50, marker='o')
    plt.plot(x_fit, y_fit, label='Fitted Parabola', color='red', linewidth=2)
    plt.axvline(x=b_fit, color='grey', linestyle='--', label=f'x = {b_fit:.4f}')
    plt.title('MZ Fit')
    plt.xlabel('mz')
    plt.ylabel('sig')
    plt.legend()
    plt.savefig("qc_%6.4f.pdf"%acc_mz)
    plt.close()

    return acc_mz

#introduction
# generate reference spectrum, ions with modification
# PEPT*IDE

#std
# IDE DE E
# PEP PE P
#mod
# PEPT*, PEPT*I, PEPT*ID
# T*IDE, PT*IDE, EPT*IDE

#usage: file scan peptide

parser = argparse.ArgumentParser(description='Analyze fragmentation patterns from mass spectrometry data')
parser.add_argument('input_file', help='Path to input file (mgf/mzML/mzXML)')
parser.add_argument('scan_number', type=int, help='Scan number to analyze')
parser.add_argument('peptide_sequence', help='Peptide sequence to match (may contain * for modifications)')
args = parser.parse_args()

full_path = args.input_file
scan = args.scan_number
peptide = args.peptide_sequence
if not peptide.replace("*", "").isalpha():
    sys.exit(f"Error: Invalid peptide sequence '{peptide}' - must contain only letters and optional '*' modification marker")

#generate ideal fragments
modpeptide = peptide.replace("*", "[+5000]") #digital labeling
prots = proforma.parse(modpeptide)
frags = fa.get_theoretical_fragments( prots[0], "aby" ) #, neutral_losses={"H2O": -18.010565}

std_mz = []
mod_mz = []
for tag, mz in frags:
    if mz > 5000:
        mod_mz.append(mz - 5000)
    else:
        std_mz.append(mz)

#generate predicted intensity with ms2pip
ms2pip_sp = SinglePrediction()
pred_mz, pred_int, _ = ms2pip_sp.predict(peptide.replace("*", ""), "-", 2, model="CID")
#match intensity with mod
max_i = np.max( pred_int )
min_i = np.min( pred_int )
scale = 2.0 #
tol = 1e-6
std_int = []
for mz in std_mz:
    d2 = [(m-mz)**2 for m in pred_mz]
    ndx = np.argmin(d2)
    min_d2 = d2[ndx]
    if min_d2 < tol:
        hit_int = pred_int[ndx]
    else:
        #should be a-ion
        hit_int = max_i/10.0
    if hit_int < max_i/scale:
        hit_int = max_i/scale
    std_int.append(hit_int)
mod_int = []
for mz in mod_mz:
    d2 = [(m-mz)**2 for m in pred_mz]
    ndx = np.argmin(d2)
    min_d2 = d2[ndx]
    if min_d2 < tol:
        hit_int = pred_int[ndx]
    else:
        #should be a-ion
        hit_int = max_i/10.0
    if hit_int < max_i/scale:
        hit_int = max_i/scale
    mod_int.append(hit_int)

val_std = np.zeros([nbin+1])
for center, h in zip(std_mz, std_int):
    val_std += h/max_i * gaussian(xbin, alpha, center)
val_mod = np.zeros([nbin+1])
for center, h in zip(mod_mz, mod_int):
    val_mod += h/max_i * gaussian(xbin, alpha, center)

#read experiment MS2
file_name = full_path.split('/')[-1]
file_type = file_name.split('.')[-1]
file_name = file_name[:-len(file_type)-1]
try:
    if file_type == "mgf":
        spectrum = get_spectrum_from_mgf(full_path, scan)
    elif file_type == "mzML":
        spectrum = get_spectrum_from_mzML(full_path, scan)
    elif file_type == "mzXML":
        spectrum = get_spectrum_from_mzXML(full_path, scan)
    else:
        sys.exit(f"Error: Unsupported file type '{file_type}'. Supported types: mgf, mzML, mzXML")
except Exception as e:
    sys.exit(f"Error reading spectrum from {full_path} (scan {scan}): {str(e)}")

if not hasattr(spectrum, 'mz') or not hasattr(spectrum, 'intensity'):
    sys.exit(f"Error: Invalid spectrum data in {full_path} (scan {scan}) - missing mz or intensity arrays")

#conv plot for exp MS2
mz = spectrum.mz
#max_it = np.max(spectrum.intensity)
sorted_it = sorted(spectrum.intensity, reverse=True)
top1 = sorted_it[0]
top2 = sorted_it[1]
top3 = sorted_it[2]
max_it = top3 #normlize using 3rd peak
raw_it = spectrum.intensity / max_it
#cap
it = [1.0 if i>1 else i for i in raw_it]
val = np.zeros([nbin+1])
for center, h in zip(mz, it):
    val += h * gaussian(xbin, alpha, center)

#first round, locate standard ions
sig_std = signal.correlate(val, val_std, mode='full')
sig_std = sig_std[halfL:-halfL][::-1]
x_i = np.argmax(sig_std)
if x_i != halfL:
    print(x_i-halfL)
    print("Warning! check the peptide, it does not match with the MS2")
else:
    print("exp_ion_number_raw:", len(mz))
    print("Peptide fit!")
    print("best_corr_score:", sig_std[x_i])

#remove std ions
mz_rm_std = remove_ions_from( mz, std_mz )
it_rm_std = extract_intensity_from( mz, it, mz_rm_std )
val = np.zeros([nbin+1])
for center, h in zip(mz_rm_std, it_rm_std):
    val += h * gaussian(xbin, alpha, center)
#check second round
sig_std = signal.correlate(val, val_std, mode='full')
sig_std = sig_std[halfL:-halfL][::-1]
x_i = np.argmax(sig_std)
h_std = sig_std[x_i]
x_p = xbin[x_i] - x_mid
print("std_global_shift:", "(%4.2f, %4.2f)"%(-x_p, h_std))
print("second_corr_score:", sig_std[x_i])
##top10
#shift = []
#for i in range(10):
#    x_i = np.argmax(sig_std)
#    x_p = xbin[x_i] - x_mid
#    shift.append(x_p)
#    h_p0 = sig_std[x_i]
#    sig_std[x_i] = 0.0
#    print("top", i+1, "(%4.2f, %4.2f)"%(-x_p, h_p0))

##check mod ions
sig_mod = signal.correlate(val, val_mod, mode='full')
sig_mod = sig_mod[halfL:-halfL][::-1]
##global
x_i = np.argmax(sig_mod)
h_mod = sig_mod[x_i]
h_std = sig_std[x_i]
x_p = xbin[x_i] - x_mid
x_top1 = x_p
print("exp_ion_number_1:", len(mz_rm_std))
print("mod_global_shift_1", "(%4.2f, %4.2f, %4.2f)"%(-x_p, h_mod, h_std))
##top10
#shift = []
#for i in range(10):
#    x_i = np.argmax(sig_mod[:halfL])
#    x_p = xbin[x_i] - x_mid
#    shift.append(x_p)
#    h_p = sig_mod[x_i]
#    h_p0 = sig_std[x_i]
#    sig_mod[x_i] = 0.0
#    print("top", i+1, "(%4.2f, %4.2f, %4.2f)"%(-x_p, h_p, h_p0))
#fit accurate top 1 mz
x_i = np.argmax(sig_mod[:halfL])
x_p = xbin[x_i] - x_mid
acc_mz = x_p + fit_mz_by_sig_peak( sig_mod, x_i )
print("fit_mod_1:", -acc_mz )
final_shift = -acc_mz
if final_shift<0:
    mod_str = "[" + "%6.4f"%final_shift + "]"
else:
    mod_str = "[+" + "%6.4f"%final_shift + "]"
#print(mod_str)
modpeptide = peptide.replace("*", mod_str) #digital labeling
print( "output", file_name, file_type, scan, modpeptide )
draw_spect_pep_pdf( spectrum, modpeptide, "spect-"+file_name+"_"+str(scan)+"_mod_1.pdf" )

#remove mod ions
mz_rm_mod = remove_ions_from( mz_rm_std, mod_mz-x_top1 )
it_rm_mod = extract_intensity_from( mz_rm_std, it_rm_std, mz_rm_mod )
val = np.zeros([nbin+1])
for center, h in zip(mz_rm_mod, it_rm_mod):
    val += h * gaussian(xbin, alpha, center)
##check mod ions, again
sig_mod = signal.correlate(val, val_mod, mode='full')
sig_mod = sig_mod[halfL:-halfL][::-1]
##global
x_i = np.argmax(sig_mod)
h_mod = sig_mod[x_i]
h_std = sig_std[x_i]
x_p = xbin[x_i] - x_mid
print("exp_ion_number_2", len(mz_rm_mod))
print("mod_global_shift_2", "(%4.2f, %4.2f, %4.2f)"%(-x_p, h_mod, h_std))
##top10
#shift = []
#for i in range(10):
#    x_i = np.argmax(sig_mod[:halfL])
#    x_p = xbin[x_i] - x_mid
#    shift.append(x_p)
#    h_p = sig_mod[x_i]
#    h_p0 = sig_std[x_i]
#    sig_mod[x_i] = 0.0
#    print("top", i+1, "(%4.2f, %4.2f, %4.2f)"%(-x_p, h_p, h_p0))
#fit accurate top 1 mz
x_i = np.argmax(sig_mod[:halfL])
x_p = xbin[x_i] - x_mid
acc_mz = x_p + fit_mz_by_sig_peak( sig_mod, x_i )
print("fit_mod_2:", -acc_mz )
final_shift = -acc_mz
if final_shift<0:
    mod_str = "[" + "%6.4f"%final_shift + "]"
else:
    mod_str = "[+" + "%6.4f"%final_shift + "]"
#print(mod_str)
modpeptide = peptide.replace("*", mod_str) #digital labeling
print( "output", file_name, file_type, scan, modpeptide )
draw_spect_pep_pdf( spectrum, modpeptide, "spect-"+file_name+"_"+str(scan)+"_mod_2.pdf" )
print( "Done!" )

