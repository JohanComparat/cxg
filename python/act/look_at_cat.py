import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
import astropy.io.fits as fits
import astropy.units as u
from glob import glob
import os
from scipy.interpolate import interp1d

from astropy.cosmology import FlatLambdaCDM, Planck18
import astropy.units as u
cosmo = Planck18
V070=cosmo.comoving_volume(0.7) *9_000*np.pi/129600./1e9
V095=cosmo.comoving_volume(0.95)*9_000*np.pi/129600./1e9
V109=cosmo.comoving_volume(1.09)*9_000*np.pi/129600./1e9
V135=cosmo.comoving_volume(1.35)*9_000*np.pi/129600./1e9

eu_dir = '/home/comparat/sf_Shared/data/Euclid'
hdu = fits.open( os.path.join( eu_dir, 'EUC_LE3_VMPZ-ID_HPDEPTHMAP-CONCAT-20241112T111230.257510Z_0.11.fits'))

area_pix = hp.nside2pixarea(hdu[1].header['NSIDE'], degrees = True)

depth_5s = hdu[1].data['WEIGHT'] - 2.5*np.log10(2)
depth_5s.sort()
clean_depth = depth_5s[depth_5s>0]
areas = np.ones_like(clean_depth) * area_pix

a0,a1 = hp.pix2ang(hdu[1].header['NSIDE'], hdu[1].data['PIXEL'][depth_5s>0], lonlat=True, nest=True)

hd2 = fits.open( os.path.join( eu_dir, 'ezgal.euc_VIS.speclite.bc03_exp_0.1_z_0.02_salp.model.mstar.cosmo.esutil.z3.fit'))
m1,z1 = np.loadtxt(os.path.join( eu_dir, 'zvlim-magVIS_0p2Lstar.ascii'), unpack = True)
itp02 = interp1d(np.hstack((m1, 30)), np.hstack((z1,2)))
itp02_r = interp1d(np.hstack((z1,2)), np.hstack((m1, 30)))
m1_04,z1_04 = np.loadtxt(os.path.join( eu_dir, 'zvlim-magVIS_0p4Lstar.ascii'), unpack = True)
itp04 = interp1d(np.hstack((m1_04, 30)), np.hstack((z1_04,2)))
itp04_r = interp1d(np.hstack((z1_04,2)), np.hstack((m1_04, 30)))

zvlim02 = itp02(clean_depth)
zvlim04 = itp04(clean_depth)

areas = np.cumsum(areas)[::-1]

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 12})
import matplotlib.pyplot as plt

s0 = (a0<70)
plt.scatter(a0[s0][::3], a1[s0][::3], c=clean_depth[s0][::3], s=1, rasterized=True, cmap='autumn')
plt.colorbar(label=r'VIS mag depth $5\sigma-0.75$')
plt.xlabel('lon healpy [deg]')
plt.ylabel('lat healpy [deg]')
#plt.yscale('log')
plt.xlim((49,69))
plt.ylim((-55,-22))
plt.savefig('lon-lat-depth-0.png')
plt.clf()

s0 = (a0>70)
plt.scatter(a0[s0], a1[s0], c=clean_depth[s0], s=1, rasterized=True, cmap='autumn')
plt.colorbar(label=r'VIS mag depth $5\sigma-0.75$')
plt.xlabel('lon healpy [deg]')
plt.ylabel('lat healpy [deg]')
#plt.yscale('log')
plt.xlim((260,280))
plt.ylim((60,72))
plt.savefig('lon-lat-depth-1.png')
plt.clf()

plt.plot(zvlim02, areas, label='0.2L*', color='b')
plt.plot(zvlim04, areas, label='0.4L*', color='r')
plt.axvline(itp02(25), ls='dotted', color='b', label='mag=24, 24.5, 25')
plt.axvline(itp02(24.5), ls='dotted', color='b')
plt.axvline(itp02(24), ls='dotted', color='b')
plt.axvline(itp04(25), ls='dotted', color='r')
plt.axvline(itp04(24.5), ls='dotted', color='r')
plt.axvline(itp04(24), ls='dotted', color='r')
plt.legend(loc=3)
plt.xlabel(r'z$_{vlim}$')
plt.ylabel('area [sq. deg.]')
#plt.yscale('log')
plt.ylim((0,70))
plt.savefig('vlim-area.png')
plt.clf()

f_scale = 10000/areas[0]
plt.plot(zvlim02, f_scale*areas, label='0.2L*', color='b')
plt.plot(zvlim04, f_scale*areas, label='0.4L*', color='r')
plt.axvline(itp02(25), ls='dotted', color='b', label='mag=24, 24.5, 25')
plt.axvline(itp02(24.5), ls='dotted', color='b')
plt.axvline(itp02(24), ls='dotted', color='b')
plt.axvline(itp04(25), ls='dotted', color='r')
plt.axvline(itp04(24.5), ls='dotted', color='r')
plt.axvline(itp04(24), ls='dotted', color='r')
plt.legend(loc=3)
plt.xlabel(r'z$_{vlim}$')
plt.ylabel('area [sq. deg.], Q1 scaled to 10kdeg2')
#plt.yscale('log')
plt.ylim((0,11_000))
plt.tight_layout()
plt.savefig('vlim-area-scaled10k.png')
plt.clf()

plt.plot(zvlim02, areas/areas.max(), label='0.2L*', color='b')
plt.plot(zvlim04, areas/areas.max(), label='0.4L*', color='r')
i2 = interp1d(areas/areas.max(), zvlim02)
i4 = interp1d(areas/areas.max(), zvlim04)
plt.axhline(0.9, ls='dotted', color='k')
plt.axvline(i2(0.9), ls='dotted', color='b', label='z='+str(np.round(i2(0.9),2))+', mag='+str(np.round(itp02_r(i2(0.9)),2)) )
plt.axvline(i4(0.9), ls='dotted', color='r', label='z='+str(np.round(i4(0.9),2))+', mag='+str(np.round(itp04_r(i4(0.9)),2)) )
#plt.axvline(itp02(25), ls='dotted', color='b', label='mag=24, 24.5, 25')
#plt.axvline(itp02(24.5), ls='dotted', color='b')
#plt.axvline(itp02(24), ls='dotted', color='b')
#plt.axvline(itp04(25), ls='dotted', color='r')
#plt.axvline(itp04(24.5), ls='dotted', color='r')
#plt.axvline(itp04(24), ls='dotted', color='r')
plt.legend(loc=3)
plt.xlabel(r'z$_{vlim}$')
plt.ylabel('area fraction')
#plt.yscale('log')
plt.ylim((0.3,1.02))
plt.savefig('vlim-area-percent.png')
plt.clf()

plt.plot(m1,z1, label='0.2L*', color='b')
plt.plot(m1_04,z1_04, label='0.4L*', color='r')
plt.axvline(24, label='mag=24, 25', ls='dashed', color='k')
plt.axvline(24.75, label='mag=24, 24.75, 25', ls='dashed', color='k')
#plt.axvline(25, ls='dashed', color='k')
plt.axhline(1.35, ls='dotted', color='k')
plt.axhline(1.1, ls='dotted', color='k')
#plt.axhline(itp02(25), ls='dashed', color='b')
#plt.axhline(itp02(24.5), ls='dashed', color='b')
#plt.axhline(itp02(24), ls='dashed', color='b')
#plt.axhline(itp04(25), ls='dashed', color='r')
#plt.axhline(itp04(24.5), ls='dashed', color='r')
#plt.axhline(itp04(24), ls='dashed', color='r')
plt.hist(clean_depth,bins=np.arange(20, 26, 0.1), weights=np.ones_like(clean_depth)*3/len(clean_depth))
plt.hist(clean_depth,bins=np.arange(20, 26, 0.1), weights=np.ones_like(clean_depth)/len(clean_depth), cumulative=True, histtype='step', lw=2, label='Euclid Q1 CDF')
plt.xlim((16,27))
plt.ylim((0.0,1.8))
plt.grid()
plt.legend(loc=2)
plt.ylabel(r'z$_{vlim}$')
plt.xlabel(r'Euclid VIS mag lim $5\sigma-0.75$')
#plt.yscale('log')
plt.savefig('vlim-magVIS.png')
plt.clf()
