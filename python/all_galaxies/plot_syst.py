print('starts')
import time
t0 = time.time()
import os, sys, glob
#import numpy as n
import numpy as np
from astropy.io import fits
import healpy
#import healsparse
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

#import healsparse as hsp
import healpy as hp
import pandas as pd
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord

nl = lambda selection : len(selection.nonzero()[0])
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT

NSIDE_val = int(sys.argv[1])
topdir    = sys.argv[2]
p2_data   = sys.argv[3]
p2_random = sys.argv[4]
p_2_2PCF  = os.path.join(topdir, sys.argv[5])
basename = sys.argv[5].split('-')[0][9:]
z_min = float(sys.argv[6])
z_max = float(sys.argv[7])

NSIDE = NSIDE_val

syst_dir  =   os.path.join(topdir, basename+'_syst')
base_out = os.path.join( syst_dir, basename+'_NSIDE_'+str(NSIDE)+'_sys' )
base_out_R = os.path.join( syst_dir, basename+'_NSIDE_'+str(NSIDE)+'_sys_RAND' )

area_pix = healpy.nside2pixarea(NSIDE, degrees = True)

fig_dir  =  os.path.join( '../figures', basename+'_syst', 'NSIDE_'+str(NSIDE).zfill(4) )
os.system('mkdir -p '+ fig_dir)

deltaMag = 0.05
deltaEBV = 0.002
deltaPSF = 0.05
deltaSD = 0.05

print(base_out)
p2_mean = base_out + '.mean.csv'
p2_sum  = base_out + '.sum.csv'
df_mean = pd.read_csv(p2_mean)
df_sum = pd.read_csv(p2_sum)

print(base_out_R)
p2_mean = base_out_R + '.mean.csv'
p2_sum  = base_out_R + '.sum.csv'
df_mean_R = pd.read_csv(p2_mean)
df_sum_R = pd.read_csv(p2_sum)

def plot_syst( sys_key='log10_stellar_density', NSIDE=NSIDE, name='',
               title='',dlognh = 0.1, y0=-0.25, y1=0.25,
               df_mean=df_mean, df_sum=df_sum,
               df_mean_R=df_mean_R, df_sum_R=df_sum_R):
    #
    # FIGURE EBV
    #
    print(sys_key)
    fig_name = os.path.join( fig_dir, basename + '_NSIDE_'+str(NSIDE)+'_syst_' + sys_key + name +'.png')

    mean_density =  np.mean(df_sum_R['w_i']/area_pix)
    p.figure(0, (7,5) )
    #
    # running mean in bins by hand
    # need to write the proper pandas command so it is faster
    #
    #boundaries = 10**np.arange( np.log10(df_mean_R[sys_key].min()), np.log10(df_mean_R[sys_key].max())+dlognh, dlognh )
    V0 = np.min(df_mean_R[sys_key][np.isinf(df_mean_R[sys_key])==False])
    V1 = np.max(df_mean_R[sys_key][np.isinf(df_mean_R[sys_key])==False])
    boundaries = np.arange( V0, V1 + dlognh, dlognh )
    x_low, x_mid, x_med, x_hig = [], [], [], []
    y_mid, y_std = [], []
    for ii in np.arange(len(boundaries)-1):
        s_nh = (df_mean_R[sys_key] >= boundaries[ii]) & ( df_mean_R[sys_key] < boundaries[ii+1])
        x_low .append( boundaries[ii]   )
        x_hig .append( boundaries[ii+1] )
        x_mid .append( 0.5 * ( boundaries[ii+1] + boundaries[ii] ) )
        x_med .append( np.mean(df_mean_R[sys_key][s_nh]) )
        y_mid .append( np.sum(df_sum_R['w_i'][s_nh])/len(df_sum_R['w_i'][s_nh]) /area_pix / mean_density -1 )
        std_frac = np.std(df_sum_R['w_i'][s_nh])/np.sum(df_sum_R['w_i'][s_nh])
        y_std .append( ( 1/np.sum(df_sum_R['w_i'][s_nh]) + std_frac**2 ) **0.5 )

    x_low = np.array( x_low )
    x_hig = np.array( x_hig )
    x_mid = np.array( x_mid )
    x_med = np.array( x_med )
    y_mid = np.array( y_mid )
    y_std = np.array( y_std )
    #print("y_mid", y_mid)
    #print("y_std", y_std)
    y_up = abs( ( 1 + y_std ) * y_mid )
    y_lo = abs( ( 1 - y_std ) * y_mid )
    ok = (y_mid>-1000)
    #print(y_up[ok] - y_mid[ok], y_mid[ok] - y_lo[ok])
    #p.errorbar(x_mid, y_mid, yerr = y_std * y_mid, xerr=[x_hig-x_mid, x_mid-x_low], lw=2, fmt='', ls='')
    p.errorbar(x_mid[ok], y_mid[ok], yerr = abs(y_std[ok]*y_mid[ok]), xerr=[x_hig[ok]-x_mid[ok], x_mid[ok]-x_low[ok]],
               lw=3, marker='*', fmt='', ls='', label='Randoms')
    RR_interp = interp1d( np.copy(x_mid), np.copy(y_mid) )
    RR_std_interp = interp1d( np.copy(x_mid), np.copy(y_std)/np.copy(y_mid) )

    mean_density =  np.mean(df_sum['w_i']/area_pix)
    #
    # running mean in bins by hand
    # need to write the proper pandas command so it is faster
    #
    #boundaries = 10**np.arange( np.log10(df_mean[sys_key].min()), np.log10(df_mean[sys_key].max())+dlognh, dlognh )
    V0 = np.min(df_mean[sys_key][np.isinf(df_mean[sys_key])==False])
    V1 = np.max(df_mean[sys_key][np.isinf(df_mean[sys_key])==False])
    boundaries = np.arange( V0, V1 + dlognh, dlognh )
    x_low, x_mid, x_med, x_hig = [], [], [], []
    y_mid, y_std = [], []
    for ii in np.arange(len(boundaries)-1):
        s_nh = (df_mean[sys_key] >= boundaries[ii]) & ( df_mean[sys_key] < boundaries[ii+1])
        x_low .append( boundaries[ii]   )
        x_hig .append( boundaries[ii+1] )
        x_mid .append( 0.5 * ( boundaries[ii+1] + boundaries[ii] ) )
        x_med .append( np.mean(df_mean[sys_key][s_nh]) )
        y_mid .append( np.sum(df_sum['w_i'][s_nh])/len(df_sum['w_i'][s_nh]) /area_pix / mean_density -1 )
        std_frac = np.std(df_sum['w_i'][s_nh])/np.sum(df_sum['w_i'][s_nh])
        y_std .append( ( 1/np.sum(df_sum['w_i'][s_nh]) + std_frac**2 ) **0.5 )

    x_low = np.array( x_low )
    x_hig = np.array( x_hig )
    x_mid = np.array( x_mid )
    x_med = np.array( x_med )
    y_mid = np.array( y_mid )
    y_std = np.array( y_std )
    #print("y_mid", y_mid)
    #print("y_std", y_std)
    y_up = abs( ( 1 + y_std ) * y_mid )
    y_lo = abs( ( 1 - y_std ) * y_mid )
    ok = (y_mid>-1000)
    #print(y_up[ok] - y_mid[ok], y_mid[ok] - y_lo[ok])
    #p.errorbar(x_mid, y_mid, yerr = y_std * y_mid, xerr=[x_hig-x_mid, x_mid-x_low], lw=2, fmt='', ls='')
    p.errorbar(x_mid[ok], y_mid[ok], yerr = abs(y_std[ok]*y_mid[ok]), xerr=[x_hig[ok]-x_mid[ok], x_mid[ok]-x_low[ok]],
               lw=3, marker='*', fmt='', ls='', label='Galaxies')
    DD_interp = interp1d(np.copy(x_mid), np.copy(y_mid))
    DD_std_interp = interp1d(np.copy(x_mid), np.copy(y_std)/np.copy(y_mid))
    out_H = np.histogram(df_mean[sys_key], bins=boundaries)#, density = True)
    y_hist = out_H[0]/np.sum(out_H[0])
    y_hist_cumul = np.cumsum(out_H[0])/np.sum(out_H[0])
    x_hist = ( out_H[1][1:] + out_H[1][:-1] ) /2.
    cms = interp1d(y_hist_cumul, x_hist)
    p.step(x_hist, 0.15 * y_hist/y_hist.max(), where='mid', label='PDF' )
    p.axvline(cms(0.1), ls='dashed', color='r', label='10-90%' )
    # p.axvline(cms(0.99), ls='dashed', color='r' )
    # p.axvline(cms(0.05), ls='dotted', color='m', label='5-95%' )
    p.axvline(cms(0.9), ls='dashed', color='r' )
    print(cms(0.1),cms(0.9))

    x_min = np.max([DD_interp.x.min(),RR_interp.x.min()])
    x_max = np.min([DD_interp.x.max(),RR_interp.x.max()])
    x_array_diff = np.arange(x_min, x_max, (x_max-x_min)/100 )
    y_diff = DD_interp(x_array_diff)-RR_interp(x_array_diff)
    y_diff_err = y_diff*(DD_std_interp(x_array_diff)**2-RR_std_interp(x_array_diff)**2)**0.5
    p.errorbar(x_array_diff, y_diff, yerr=y_diff_err, fmt='o', lw=3, color='k', label=r'$\Delta$ data-random')
    p.xlabel(sys_key)
    p.ylabel(r"$N/\bar{N}-1$ ")
    #if sys_key=='EBV':
        #p.xlim((0.0, 0.6))
        #p.axvline(0.15, ls='dashed')
    p.ylim((y0, y1))
    p.xlim((V0, V1))
    p.grid()
    p.legend()
    p.title(title+', NSIDE=' + str(NSIDE))
    #p.xscale('log')
    p.tight_layout()
    p.savefig( fig_name )
    p.clf()
    print(fig_name)

# stellar density


deltaMag = 0.05
deltaEBV = 0.002
deltaPSF = 0.05
deltaSD = 0.05

z_str = '0.1<z<0.2'

plot_syst( 'log10_stellar_density',
           title=z_str,
           dlognh = deltaSD,
            df_mean=df_mean, df_sum=df_sum,
            df_mean_R=df_mean_R, df_sum_R=df_sum_R)#, y0=-1, y1=1 )

