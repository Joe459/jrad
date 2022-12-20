use_MPI = True
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import csv
import sys
if use_MPI:
	from mpi4py.futures import MPIPoolExecutor
	from mpi4py import MPI
import time
from scipy.signal import find_peaks, peak_widths, peak_prominences, savgol_filter, medfilt
from scipy.optimize import curve_fit
from functools import partial
import warnings


def profile(file_name):
	
	dat = np.genfromtxt('./bunch-prof/' + file_name)
	print(dat)

if __name__ == '__main__':
	
	if use_MPI:			
		with MPIPoolExecutor() as ex:
			ex.map(profile,os.listdir('./bunch-prof/'))
	else:
		[profile(file_name) for file_name in os.listdir('./bunch-prof/')]	
