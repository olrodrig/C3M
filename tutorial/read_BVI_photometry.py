import numpy as np

def read_BVI_photometry(input_file):
  t  = np.genfromtxt(input_file, usecols=0)
  B  = np.genfromtxt(input_file, usecols=1)
  eB = np.genfromtxt(input_file, usecols=2)
  V  = np.genfromtxt(input_file, usecols=3)
  eV = np.genfromtxt(input_file, usecols=4)
  I  = np.genfromtxt(input_file, usecols=5)
  eI = np.genfromtxt(input_file, usecols=6)
  
  photometry = {'t': t, 'B':B, 'eB':eB, 'V':V, 'eV':eV, 'I':I, 'eI':eI}
  
  return photometry
