"""
python compare_wp.py 0 40
python compare_wp.py 1 40
python compare_wp.py 2 40
python compare_wp.py 3 40
python compare_wp.py 4 40
python compare_wp.py 5 40
python compare_wp.py 6 40
python compare_wp.py 7 40

python compare_wp.py 0 60
python compare_wp.py 1 60
python compare_wp.py 2 60
python compare_wp.py 3 60
python compare_wp.py 4 60
python compare_wp.py 5 60
python compare_wp.py 6 60
python compare_wp.py 7 60

python compare_wp.py 0 100
python compare_wp.py 1 100
python compare_wp.py 2 100
python compare_wp.py 3 100
python compare_wp.py 4 100
python compare_wp.py 5 100
python compare_wp.py 6 100
python compare_wp.py 7 100

"""
print('+'*100)
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import os, sys
import glob
import numpy as np
from astropy.table import Table, Column, vstack, hstack
from scipy.optimize import curve_fit
import healpy
from scipy.interpolate import interp1d
colors = ["#67E568","#FFF000","#FFB62B","#E56124",
		  "#E53E30",
		  "#7F2353","#F911FF","#9F8CA6","#257F27","#08420D"]

import astropy.units as u
import astropy.constants as cc
from astropy.cosmology import FlatLambdaCDM
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT


fig_dir  ='../../figures/'
dat_dir  ='../../data/'
deg_to_rad = np.pi/180.

cl_sample = sys.argv[1] # 'S0'
p2_all_dat  = np.array(glob.glob(os.path.join(dat_dir, 'Counts_DATA_'+cl_sample+'_RAND_gr_*_z_*.RSdata.npy')))

all_z         = []
all_colors    = []
all_sn_1Mpc   = []
all_sn_100kpc = []

color_step = 0.05

RES = {}
RES_MAX = []
RES_MAX_100kpc = []
for jj, p2_data in enumerate(p2_all_dat):
    RES[jj] = np.load(p2_data, allow_pickle='TRUE').item()
    all_z         .append(RES[jj]['z_bar'] * np.ones_like(RES[jj]['color_grid']))
    all_colors    .append(RES[jj]['color_grid']+color_step/2.)
    all_sn_1Mpc   .append(RES[jj]['sumsn_1Mpc_all'])
    all_sn_100kpc .append(RES[jj]['sumsn_100kpc_all'])
    i_max = np.argmax(RES[jj]['sumsn_1Mpc_all'])
    RES_MAX.append([RES[jj]['z_bar'], RES[jj]['color_grid'][i_max]+color_step/2.])
    i_max = np.argmax(RES[jj]['sumsn_100kpc_all'])
    RES_MAX_100kpc.append([RES[jj]['z_bar'], RES[jj]['color_grid'][i_max]+color_step/2.])

RES_MAX = np.transpose(RES_MAX)
RES_MAX_100kpc = np.transpose(RES_MAX_100kpc)

all_z         = np.hstack((all_z        ))
all_colors    = np.hstack((all_colors   ))
all_sn_1Mpc   = np.hstack((all_sn_1Mpc  ))
all_sn_100kpc = np.hstack((all_sn_100kpc))
all_sn_1Mpc[all_sn_1Mpc<=0]=0.001
all_sn_100kpc[all_sn_100kpc<=0]=0.001
RS = {}
RS['all_sn_100kpc'] = all_sn_100kpc
RS['all_sn_1Mpc']   = all_sn_1Mpc
RS['all_colors']    = all_colors
RS['all_z']         = all_z

RS_model = Table.read( os.path.join( os.environ['GIT_STMOD_DATA'], 'data', 'models', 'model_GAL', 'legacy_dr10_south_v0.3_grz_z_cal_zspec_redgals_model.fit') )
z_RS = np.hstack(( RS_model['nodes'][0] ))
gr_RS = np.hstack(( np.array(RS_model['meancol']).T[0].T[0] ))
gr_RS_sigma = np.hstack(( np.array(RS_model['meancol_scatter']).T[0].T[0] ))

# plot CT

p2_fig  = os.path.join(fig_dir, 'RS_gr_z_SN100kpc_'+cl_sample+'.png')
plt.figure(11, (7.4,7.5))
plt.scatter(RS['all_z'], RS['all_colors'], s=100, marker='s', c=np.log10(RS['all_sn_100kpc']) , cmap=pl.cm.cool_r, vmin=0, label='GalxClusters', alpha=0.3)#, vmax=15)
plt.plot(z_RS, gr_RS, 'k--', lw=2, label='model Kluge 24')
plt.plot(z_RS, gr_RS-gr_RS_sigma, 'k:', lw=2)
plt.plot(z_RS, gr_RS+gr_RS_sigma, 'k:', lw=2)
plt.plot(RES_MAX[0], RES_MAX[1], ls='--',lw=2,label='highest S/N 1Mpc',c='darkred')
plt.plot(RES_MAX_100kpc[0], RES_MAX_100kpc[1], ls='--',lw=2,label='highest S/N 100kpc',c='darkgreen')

plt.colorbar(label=r'log10 $\Sigma [S/N(r<100kpc)]$')
plt.xlim((RS['all_z'].min()-0.01, RS['all_z'].max()+0.01))
plt.ylim((RS['all_colors'].min()-0.01, RS['all_colors'].max()+0.01))
plt.ylabel(r'color $g-r$')
plt.xlabel(r"$z$")
plt.title('clusters '+cl_sample )
plt.legend(loc=4, fontsize=12, framealpha=0.95)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

p2_fig  = os.path.join(fig_dir, 'RS_gr_z_SN1Mpc_'+cl_sample+'.png')
plt.figure(11, (7.4,7.5))
plt.scatter(RS['all_z'], RS['all_colors'], s=100, marker='s', c=np.log10(RS['all_sn_1Mpc']) , cmap=pl.cm.cool_r, vmin=0, label='GalxClusters', alpha=0.3)#, vmax=15)
plt.plot(z_RS, gr_RS, 'k--',lw=2,label='model Kluge 24')
plt.plot(z_RS, gr_RS-gr_RS_sigma, 'k:', lw=2)
plt.plot(z_RS, gr_RS+gr_RS_sigma, 'k:', lw=2)
plt.plot(RES_MAX[0], RES_MAX[1], ls='--',lw=2,label='highest S/N 1Mpc',c='darkred')
plt.plot(RES_MAX_100kpc[0], RES_MAX_100kpc[1], ls='--',lw=2,label='highest S/N 100kpc',c='darkgreen')
plt.colorbar(label=r'log10 $\Sigma [S/N(r<1Mpc)]$')
plt.xlim((RS['all_z'].min()-0.01, RS['all_z'].max()+0.01))
plt.ylim((RS['all_colors'].min()-0.01, RS['all_colors'].max()+0.01))
plt.ylabel(r'color $g-r$')
plt.xlabel(r"$z$")
plt.title('clusters '+cl_sample )
plt.legend(loc=4, fontsize=12, framealpha=0.95)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)




