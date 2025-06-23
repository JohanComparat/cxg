import time
t0 = time.time()
import numpy as np
import numpy as n
import healpy
from astropy.table import Table, Column, vstack, hstack
import sys, os, glob

import astropy.units as u
import astropy.constants as cc
import astropy.io.fits as fits

speed_light = cc.c.to(u.km/u.s).value
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
import Corrfunc
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from Corrfunc.io import read_catalog
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.theory.DD import DD
from numpy.random import choice

print(sys.argv)
print(len(sys.argv))
pimax = int(sys.argv[1])
topdir    = sys.argv[2]
p2_data   = sys.argv[3]
p2_random = sys.argv[4]
p_2_2PCF  = os.path.join(topdir, sys.argv[5])
basename = sys.argv[5].split('-')[0][9:]
z_min = float(sys.argv[6])
z_max = float(sys.argv[7])
NSIDE_val = int(sys.argv[8])
out_subfolder = sys.argv[9]
bin_size = float(sys.argv[10])

JK_dir = os.path.join(topdir, out_subfolder)
os.system('mkdir -p '+JK_dir)
print(JK_dir)

def tabulate_wprp_clustering_noW(RA, DEC, Z, rand_RA , rand_DEC, rand_Z, out_file='test.fits', pimax = 100.0, bin_size=bin_size ):
	CZ = Z * speed_light
	rand_CZ = rand_Z * speed_light
	N = len(RA)
	rand_N = len(rand_RA)
	#print(N, rand_N, out_file, time.time()-t0)
	bins = 10**np.arange(-2.0, 1.81, bin_size)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 32
	# Auto pairs counts in DD
	autocorr=1
	DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
								np.array(list(RA)).astype('float'),
								np.array(list(DEC)).astype('float'),
								np.array(list(CZ)).astype('float'))
	print('DD', DD_counts['npairs'])
	# Auto pairs counts in DR
	autocorr=0
	DR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							np.array(list(RA)).astype('float'),
							np.array(list(DEC)).astype('float'),
							np.array(list(CZ)).astype('float'),
							RA2 =rand_RA.astype('float'),
							DEC2=rand_DEC.astype('float'),
							CZ2 =rand_CZ.astype('float'))
	print('DR', DR_counts['npairs'])
	# Auto pairs counts in RR
	autocorr=1
	RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							rand_RA.astype('float'),
							rand_DEC.astype('float'),
							rand_CZ.astype('float'))
	print('RR',RR_counts['npairs'])
	# All the pair counts are done, get the angular correlation function
	wp = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts, nbins, pimax)
	t = Table()
	t.add_column(Column(data = bins[:-1], name='rp_min', unit='Mpc'  ) )
	t.add_column(Column(data = bins[1:], name='rp_max', unit='Mpc'  ) )
	t.add_column(Column(data = x, name='rp_mid', unit='Mpc'  ) )
	t.add_column(Column(data = wp, name='wprp', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * N, name='N_data', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * rand_N, name='N_random', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * pimax, name='pimax', unit=''  ) )
	#print(out_file)
	t.write(out_file, overwrite=True, format='fits')
	print(out_file, time.time()-t0, 's')

data1 = Table.read(os.path.join(topdir, p2_data ), format='fits')
rand = Table.read(os.path.join(topdir, p2_random ), format='fits')
data1=data1[(data1['Z_LAMBDA']>=z_min)&(data1['Z_LAMBDA']<=z_max)]
rand=rand[(rand['redshift']>=z_min)&(rand['redshift']<=z_max)]

NSIDE = NSIDE_val
print('='*100)
print(NSIDE)
NSIDE_str = str(NSIDE).zfill(2)
for jk_i in np.arange(1000):
	p_2_2PCF = os.path.join(JK_dir, PCF_basename + '_NSIDE_' + NSIDE_str + '_J_'+str(jk_i).zfill(4) + '.fits')
	print(jk_i, p_2_2PCF, os.path.isfile(p_2_2PCF))
	if os.path.isfile(p_2_2PCF) == False:
		U_NSIDE_list = np.unique(DDD['HPX_' + NSIDE_str])
		N_select = int(len(U_NSIDE_list)*0.9)
		pixels_keep = choice(U_NSIDE_list, size=N_select, replace=False)#, p=None)
n(RRR['HPX_' + NSIDE_str], pixels_keep)
		data1_in = np.isin(data1['HPX_' + NSIDE_str], pixels_keep)
		rand_in = np.isin(rand['HPX_' + NSIDE_str], pixels_keep)
		print(len(U_NSIDE_list), N_select, len(pixels_keep), pixels_keep[:15],
				len(data1['RA'][data1_in]), len(data1['RA'][data1_in]) / len(data1['RA']), len(rand['RA'][rand_in]) / len(rand['RA']) )
		tabulate_wprp_clustering_noW(
			data1['RA'][data1_in], data1['DEC'][data1_in], data1['Z_LAMBDA'][data1_in],
			rand['RA'][rand_in] , rand['DEC'][rand_in], rand['redshift'][rand_in],
			out_file=p_2_2PCF, pimax = 100.0)
