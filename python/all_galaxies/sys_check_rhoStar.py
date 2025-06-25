print('starts')
import time
t0 = time.time()
import os, sys, glob
import numpy as n
from astropy.io import fits
import healpy
import pandas

#from astropy_healpix import HEALPix
from astropy.coordinates import Galactic, SkyCoord
from astropy.table import Table, Column, vstack, hstack
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
import matplotlib.pylab as pl
# Planck dust map setup
#from dustmaps.planck import PlanckQuery
#from dustmaps.config import config

import healpy as hp
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord

nl = lambda selection : len(selection.nonzero()[0])
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT

import astropy.constants as cc
speed_light = cc.c.to(u.km/u.s).value

NSIDE_val = int(sys.argv[1])
topdir    = sys.argv[2]
p2_data   = sys.argv[3]
p2_random = sys.argv[4]
p_2_2PCF  = os.path.join(topdir, sys.argv[5])
basename = sys.argv[5].split('-')[0][9:]
z_min = float(sys.argv[6])
z_max = float(sys.argv[7])



data1 = Table.read(os.path.join(topdir, p2_data ), format='fits')
rand = Table.read(os.path.join(topdir, p2_random ), format='fits')

NSIDES = [4, 8, 16]
for NSIDE in NSIDES:
	NSIDE_str = str(NSIDE).zfill(2)
	data1['HPX_'+NSIDE_str] = healpy.ang2pix(NSIDE, np.pi/2. - data1['DEC']*np.pi/180. , data1['RA']*np.pi/180. , nest=True)
	rand['HPX_'+NSIDE_str]  = healpy.ang2pix(NSIDE, np.pi/2. - rand['DEC']*np.pi/180. , rand['RA']*np.pi/180. , nest=True)
	mask_dir = os.path.join(os.environ['GIT_DR10W_DATA'], 'figures/BGS_wprp/HPX_'+NSIDE_str)
	REJ_file = os.path.join( mask_dir, basename+'-wprp-pimax100-bin0p05-JK100.fits-wprp-all-bin0p05_REJoutlierClipped3Sigma.fits')
	if os.path.isfile(REJ_file) :
		t = Table.read( REJ_file )
		data1['keep_HPX_'+NSIDE_str] = (np.isin(data1['HPX_'+NSIDE_str], t['hpx_idx'], invert = True))
		rand['keep_HPX_'+NSIDE_str]  = (np.isin(rand['HPX_'+NSIDE_str], t['hpx_idx'], invert = True))
	else:
		data1['keep_HPX_'+NSIDE_str] = True
		rand['keep_HPX_'+NSIDE_str]  = True

data_keep = ( data1['keep_HPX_04'] ) & ( data1['keep_HPX_08'] ) & ( data1['keep_HPX_16'] ) & (data1['BEST_Z']>=z_min) & (data1['BEST_Z']<=z_max)
rand_keep = ( rand['keep_HPX_04'] ) & ( rand['keep_HPX_08'] ) & ( rand['keep_HPX_16'] )& (rand['Z']>=z_min) & (rand['Z']<=z_max)

t0 = time.time()
print('before masking ND, NR = ', len(data1), len(rand))
DDD = data1[ (data_keep) ]
RRR = rand[  (rand_keep) ]

N_gal = len(DDD)
N_RD = len(RRR)
RRR['Z'] = n.tile(DDD['BEST_Z'], int(N_RD*1./N_gal)+1)[:N_RD]

R_tF = n.random.random(N_RD)
N_R_F=5
sR = ( R_tF < N_gal * N_R_F / N_RD )
RRR = RRR[sR]

syst_dir  =   os.path.join(topdir, basename+'_syst')
os.system('mkdir -p '+syst_dir)

gaia_dir  =  '/data42s/comparat/firefly/gaia_cat'
t_gaia = Table.read( os.path.join( gaia_dir, 'GAIA_stellar_density_NSIDE_'+str(NSIDE).zfill(5) +'.fits' ) )
t_gaia['stellar_density'][t_gaia['stellar_density'] == 0] = 1

# EBV_MAX = 0.07
# G_GALDEPTH_MIN = 25.80
# R_GALDEPTH_MIN = 25.36
# Z_GALDEPTH_MIN = 24.40
# R_GALDEPTH_MAX = 27.22
log_rho_star_MAX = 3

base_out = os.path.join( syst_dir, basename+'_NSIDE_'+str(NSIDE)+'_sys' )
DDD['HEALPIX_IDX'] = healpy.ang2pix(NSIDE, n.pi/2. - DDD['DEC']*n.pi/180. , DDD['RA']*n.pi/180. , nest=True)
DDD['log10_stellar_density'] = np.log10(t_gaia['stellar_density'][DDD['HEALPIX_IDX']])

t2 = Table()
t2['HEALPIX_IDX'] = DDD['HEALPIX_IDX']
t2['log10_stellar_density'] = DDD['log10_stellar_density']

df = pandas.DataFrame()
df['HEALPIX_IDX'] = t2['HEALPIX_IDX']
df['log10_stellar_density'] = t2['log10_stellar_density']
df['w_i'] = 1.
df_mean     = df.groupby('HEALPIX_IDX').mean()
df_sum      = df.groupby('HEALPIX_IDX').sum()

df_mean   .to_csv(base_out +    '.mean.csv', mode='w')
df_sum    .to_csv(base_out +     '.sum.csv', mode='w')
print(base_out +    '.mean.csv')
print(base_out +    '.sum.csv')

base_out = os.path.join( syst_dir, basename+'_NSIDE_'+str(NSIDE)+'_sys_RAND' )
RRR['HEALPIX_IDX'] = healpy.ang2pix(NSIDE, n.pi/2. - RRR['DEC']*n.pi/180. , RRR['RA']*n.pi/180. , nest=True)
RRR['log10_stellar_density'] = np.log10(t_gaia['stellar_density'][RRR['HEALPIX_IDX']])

t2 = Table()
t2['HEALPIX_IDX'] = RRR['HEALPIX_IDX']
t2['log10_stellar_density'] = RRR['log10_stellar_density']

df = pandas.DataFrame()
df['HEALPIX_IDX'] = t2['HEALPIX_IDX']
df['log10_stellar_density'] = t2['log10_stellar_density']
df['w_i'] = 1.
df_mean     = df.groupby('HEALPIX_IDX').mean()
df_sum      = df.groupby('HEALPIX_IDX').sum()

df_mean   .to_csv(base_out +    '.mean.csv', mode='w')
df_sum    .to_csv(base_out +     '.sum.csv', mode='w')
print(base_out +    '.mean.csv')
print(base_out +    '.sum.csv')
