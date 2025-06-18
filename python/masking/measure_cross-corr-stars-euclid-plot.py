import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
import numpy as n
import astropy.io.fits as fits
import os, sys, glob
from os.path import join
# import pymangle as mangle
import numpy as np
import matplotlib.pyplot as p
from scipy.interpolate import interp1d

import astropy.units as u
import astropy.cosmology as cc

cosmo = cc.Planck13
from astropy.coordinates import SkyCoord
colors = ["#67E568","#FFF000","#FFB62B","#E56124","#E53E30","#7F2353","#F911FF","#9F8CA6","#257F27","#08420D"]

# from sklearn.neighbors import BallTree
from astropy.table import Table, unique
# from math import radians, cos, sin, asin, sqrt, pi

# event counts as a function of angular distance
deg_to_rad = np.pi / 180.
arcsec = 1. / 3600.
rs = np.logspace(1, 3, 40)
r_rad = deg_to_rad * rs * arcsec

# DATA
option = sys.argv[1] # "16_g_17"
fig_dir  ='../../figures/mask'
xcorr_list = np.array(glob.glob(os.path.join(fig_dir, 'xcorr_gal_' + option + '_sourceCat*_mask0.data')))
xcorr_list.sort()
max_CT = np.array([float(os.path.basename(el)[:-5].split('_')[-2]) for el in xcorr_list])
idx = np.argsort(max_CT)
xcorr_list = xcorr_list[idx][::-1]

RRcorr_list = np.array(glob.glob(os.path.join(fig_dir, 'xcorr_RR_' + option + '_sourceCat*_mask0.data')))
RRcorr_list.sort()
max_CT2 = np.array([float(os.path.basename(el)[:-5].split('_')[-2][4:]) for el in RRcorr_list])
idxR2 = np.argsort(max_CT2)
RRcorr_list = RRcorr_list[idxR2]

xcorr1_list = np.array(glob.glob(os.path.join(fig_dir, 'xcorr_gal_' + option + '_sourceCat*_mask1.data')))
xcorr1_list.sort()
max_CT1 = np.array([float(os.path.basename(el)[:-5].split('_')[-2]) for el in xcorr1_list])
idx1 = np.argsort(max_CT1)
xcorr1_list = xcorr1_list[idx1][::-1]

RRcorr1_list = np.array(glob.glob(os.path.join(fig_dir, 'xcorr_RR_' + option + '_sourceCat*_mask1.data')))
RRcorr1_list.sort()
max_CT3 = np.array([float(os.path.basename(el)[:-5].split('_')[-2][4:]) for el in RRcorr1_list])
idxR3 = np.argsort(max_CT3)
RRcorr1_list = RRcorr1_list[idxR3]

def smooth(x, window_len=11, window='hanning'):
    """
	window can be : 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
	"""
    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    if window == 'flat':
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')
    y = np.convolve(w / w.sum(), s, mode='valid')
    return y[int(window_len / 2 - 1):-int(window_len / 2)]

window_len = 10

p.figure(1, (12, 7))
for el, RRel, el1, RRel1, cc in zip(xcorr_list[::1], RRcorr_list[::1], xcorr1_list[::1], RRcorr1_list[::1], colors):
    print(el, RRel)
    x, yi = np.loadtxt(el, unpack=True)
    xR, yR = np.loadtxt(RRel, unpack=True)
    y = yi / np.mean(yR[-10:])
    y1 = smooth(y, window_len=window_len, window='blackman')
    lab = os.path.basename(el)[:-5].split('_')
    rho = np.mean(y1[-30:])
    title =  os.path.basename(el)[:-5].split('_')[-7]+'<F_VIS<' + os.path.basename(el)[:-5].split('_')[-2]
    p.plot(x, y1, label=title, lw=2, color=cc)
    el, RRel = el1, RRel1
    print(el, RRel)
    x, yi = np.loadtxt(el, unpack=True)
    xR, yR = np.loadtxt(RRel, unpack=True)
    y = yi / np.mean(yR[-10:])
    y1 = smooth(y, window_len=window_len, window='blackman')
    lab = os.path.basename(el)[:-5].split('_')
    rho = np.mean(y1[-30:])
    title =  os.path.basename(el)[:-5].split('_')[-7]+'<F_VIS<' + os.path.basename(el)[:-5].split('_')[-2]
    p.plot(x, y1, lw=2, ls='dashed', color=cc)

p.axhline(1, ls='dashed', color='k')
p.legend(frameon=True, loc=2, fontsize=10, ncol=2)
p.xlabel('separation to stars [arcseconds]')
p.ylabel('relative density of galaxies')
p.xscale('log')
# p.yscale('log')
p.title('GAIA stars, '+option)
p.savefig(os.path.join(fig_dir, 'xcorr-few-'+option+'.png'))
p.clf()


"""
cat_dir = "/home/comparat/sf_Shared/data/erosita/observations/eFEDS_c001_clean"
p_2_catalog = os.path.join(cat_dir, 'eFEDS_AGN_spec_V7.6.fits')
eRO_sources = fits.open(p_2_catalog)[2].data
mask_radius = eRO_sources['Rad']
mask_radius_ID = eRO_sources['ID']
scCt = np.log10(eRO_sources['SrcCts'])
p2_extCat = os.path.join(cat_dir, "eFEDS_c001_main_V7.3.fits")
ext_cat_i = Table.read(p2_extCat)
is_not_extended = (ext_cat_i['EXT'] == 0)  # & ( ext_cat_i['DET_LIKE']>= 8 )

p.figure(1, (8, 5))
p.axvline(9.6, color='k', ls='solid', lw=1, label='1, 3 pixels')
p.axvline(3 * 9.6, color='k', ls='solid', lw=1)
dat = []
for el, RRel in zip(xcorr_list, RRcorr_list):
    # print(el)
    x, yi = np.loadtxt(el, unpack=True)
    xR, yR = np.loadtxt(RRel, unpack=True)
    y1 = yi / np.mean(yR[-10:])
    # y1=smooth(y, window_len=window_len, window='blackman')
    lab = os.path.basename(el)[:-5].split('_')
    rho = np.mean(y1[-5:])
    below = (y1 < 1.5 * rho)
    title = lab[-3] + '-' + lab[-1]
    radius = np.min(x[below])
    print(np.min(x[below]), title)
    p.plot(x, y1 / rho, label=title, lw=2)
    # p.plot(radius, 0.8, 'k+')
    below = (y1 < 1.25 * rho)
    radius2 = np.min(x[below])
    below = (y1 < 2 * rho)
    radius3 = np.min(x[below])
    dat.append([float(lab[-3]), float(lab[-1]), radius, radius2, radius3])

p.axhline(1.25, c='grey', lw=0.7, ls='dashed', label='1.25, 1.5, 2')
p.axhline(1.5, c='grey', lw=0.7, ls='dashed')
p.axhline(2.0, c='grey', lw=0.7, ls='dashed')

p.legend(frameon=True, loc=1, fontsize=9, ncol=3)
p.xlabel('separation to source [arcseconds]')
p.ylabel('relative density of events')
p.xscale('log')
# p.xlim((0.2,10000))
# p.ylim((0.,2.))
p.title(option)
p.yscale('log')
# p.grid()
p.savefig(os.path.join(out_dir, option + 'xcorr.png'))
p.clf()

dat = np.transpose(dat)

xmin = dat[0]
xmax = dat[1]
x_data = (np.log10(xmin) + np.log10(xmax)) / 2.
mask = dat[2]
mask2 = dat[3]
mask3 = dat[4]
M0 = np.min([mask, mask2, mask3], axis=0) * 0.95
M1 = np.max([mask, mask2, mask3], axis=0) * 1.05

p.figure(2, (8, 5))
# p.plot(x_data, mask,lw=2, marker='o', label='1.1')
# p.plot(x_data, np.polyval(coeffs, x_data), ls='--')
# p.plot(x_data, mask2,lw=2, marker='o', label='1.05')
# coeffs = np.polyfit(x_data[:-3], mask2[:-3], deg=1)
# p.plot(x_data, np.polyval(coeffs, x_data), ls='--')
if option == 'EXT':
    p.plot(scCt[is_not_extended == False], mask_radius[is_not_extended == False], 'r+', label='eFEDs extended sources',
           alpha=0.3)
    coeffs = np.polyfit(scCt[(is_not_extended == False) & (scCt > 1.)],
                        mask_radius[(is_not_extended == False) & (scCt > 1.)], deg=2)
    label = r'$y=$' + str(np.round(coeffs[0], 3)) + r'$x^2+$ ' + str(np.round(coeffs[1], 3)) + r'$x+$ ' + str(
        np.round(coeffs[2], 3))
    x_poly = np.arange(0.5, 3.5, 0.01)
    p.plot(x_poly, np.polyval(coeffs, x_poly), ls='--', label=label)
    p.plot(x_poly, np.polyval(coeffs, x_poly) * 1.2, ls='--', label='x1.2')

if option == 'PS':
    p.plot(scCt[is_not_extended], mask_radius[is_not_extended], marker='+', ls='', label='eFEDs point-like sources',
           alpha=0.3)
    coeffs = np.polyfit(scCt[(is_not_extended) & (scCt > 1.)], mask_radius[(is_not_extended) & (scCt > 1.)], deg=2)
    label = r'$y=$' + str(np.round(coeffs[0], 3)) + r'$x^2+$ ' + str(np.round(coeffs[1], 3)) + r'$x+$ ' + str(
        np.round(coeffs[2], 3))
    x_poly = np.arange(0.5, 3.5, 0.01)
    p.plot(x_poly, np.polyval(coeffs, x_poly), ls='--', label=label)
    p.plot(x_poly, np.polyval(coeffs, x_poly) * 1.4, ls='--', label='x1.4')

p.errorbar(x_data, (M0 + M1) / 2., yerr=(M1 - M0) / 2., color='k', alpha=0.9, label='EVTxSRC')  # , label='1.15')
coeffs = np.polyfit(x_data, (M0 + M1) / 2., deg=2)
label = r'$y=$' + str(np.round(coeffs[0], 3)) + r'$x^2+$ ' + str(np.round(coeffs[1], 3)) + r'$x+$ ' + str(
    np.round(coeffs[2], 3))
x_poly = np.arange(0.5, 3.5, 0.01)
p.plot(x_poly, np.polyval(coeffs, x_poly), ls='-.', label=label)
# p.plot(x_data, M0, M1, lw=3, marker='o', color='m', label='xcorr 1.1')
# p.plot(x_data, mask, lw=3, marker='o', color='y', label='xcorr 1.05')#, label='1.15')
# coeffs = np.polyfit(x_data[:-3], mask3[:-3], deg=1)
# label='y='+str(np.round(coeffs[1],2))+'+'+str(np.round(coeffs[0],2))+'x'
# p.plot(x_data, np.polyval(coeffs, x_data), ls='--', label=label)
# p.plot(x_data, np.polyval(coeffs, x_data)/1.2, ls='--', label=label+'/1.5')
p.axhline(30, c='k', ls=':')
p.legend(frameon=True, loc=0, fontsize=12)
p.xlabel(r'$\log_{10}({\rm Count}\; 0.2-2.3 {\rm keV})$')
p.ylabel('masking radius [arc seconds]')
# p.xscale('log')
p.xlim((0.5, 4))
p.ylim((10., 600.))
p.yscale('log')
p.title(option)
# p.grid()
p.savefig(os.path.join(out_dir, option + 'xcorr2.png'))
p.clf()

"""