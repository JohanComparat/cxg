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

topdir    = sys.argv[1]
color_name   = sys.argv[2]
p2_data   = os.path.join(topdir, sys.argv[3])
p2_fig  = os.path.join(fig_dir, sys.argv[3][:-3]+'png')
p2_dat  = os.path.join(dat_dir, sys.argv[3][:-3]+'RSdata.npy')

#D2 = Table.read(os.path.join(C_topdir, C_p2_data ), format='fits')
#R2 = Table.read(os.path.join(C_topdir, C_p2_random ), format='fits')
#D2=D2[(D2['Z_LAMBDA']>=z_min)&(D2['Z_LAMBDA']<=z_max)]
#R2=R2[(R2['redshift']>=z_min)&(R2['redshift']<=z_max)]
#p_2_OUT='/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep/Counts_DATA_S0_RAND_gr_011_z_012.npy'
#p2_fig = os.path.join( fig_dir, 'CT_s0_gr_011z012.png')
RES = np.load(p2_data, allow_pickle='TRUE').item()
x_conversion = cosmo.kpc_proper_per_arcmin(RES['z_bar']).to(u.Mpc/u.rad).value

#RES['color_grid'] = color_grid
#RES['color_step'] = color_step
#RES['theta_grid'] = radii
#RES['theta_step'] = radius_step
#RES['DATA'] = CT_mat
#RES['RAND'] = RR_mat

sumsn_100kpc_all = []
sumsn_1Mpc_all = []
for jj, c_val in zip(np.arange(len(RES['color_grid'])), RES['color_grid']):
    diff_data = np.hstack(( RES['DATA'][jj][0], RES['DATA'][jj][1:]-RES['DATA'][jj][:-1] ))
    diff_rand = np.hstack(( RES['RAND'][jj][0], RES['RAND'][jj][1:]-RES['RAND'][jj][:-1] ))
    y_pcf = (diff_data / RES['N_D'] ) / ( diff_rand / RES['N_R'] )
    y_pcf_err = y_pcf* (1/diff_data + 1/diff_rand)**0.5
    ratio = y_pcf/y_pcf_err
    x_pcf = 0.5*( RES['theta_grid']+ np.hstack((0., RES['theta_grid'][:-1])) )
    x_pcf_up =  RES['theta_grid']
    x_pcf_lo = np.hstack((0., RES['theta_grid'][:-1]))
    sumsn_100kpc = np.sum(ratio[(~np.isnan(ratio))&(x_pcf*x_conversion<0.1)])
    sumsn_1Mpc_all.append( np.sum(ratio[(~np.isnan(ratio))&(x_pcf*x_conversion<1)]) )
    #print(sumsn_100kpc)
    sumsn_100kpc_all.append(sumsn_100kpc)

sumsn_100kpc_all = np.array(sumsn_100kpc_all)
sumsn_1Mpc_all = np.array(sumsn_1Mpc_all)
SN={}
SN['sumsn_100kpc_all'] = sumsn_100kpc_all
SN['sumsn_1Mpc_all'] = sumsn_1Mpc_all
SN['color_grid'] = RES['color_grid']
SN['z_bar'] = RES['z_bar']
np.save(p2_dat, SN)
print(p2_dat, 'written')

SN_min = 5
sel = (sumsn_100kpc_all> SN_min)
print(RES['color_grid'][sel])
if len(RES['color_grid'][sel])>0:
    #sys.exit()

    c_val_norm = (RES['color_grid'][sel] -  RES['color_grid'][sel].min() ) / ( RES['color_grid'][sel].max() - RES['color_grid'][sel].min() )
    colors = pl.cm.rainbow(c_val_norm)
    colors = pl.cm.tab20b(c_val_norm)

    # plot CT
    plt.figure(10, (7.5,5.5))

    for jj, c_val, cc in zip(np.arange(len(RES['color_grid']))[sel], RES['color_grid'][sel], colors):
        diff_data = np.hstack(( RES['DATA'][jj][0], RES['DATA'][jj][1:]-RES['DATA'][jj][:-1] ))
        diff_rand = np.hstack(( RES['RAND'][jj][0], RES['RAND'][jj][1:]-RES['RAND'][jj][:-1] ))
        y_pcf = (diff_data / RES['N_D'] ) / ( diff_rand / RES['N_R'] )
        y_pcf_err = y_pcf* (1/diff_data + 1/diff_rand)**0.5
        ratio = y_pcf/y_pcf_err
        x_pcf = 0.5*( RES['theta_grid']+ np.hstack((0., RES['theta_grid'][:-1])) )
        x_pcf_up =  RES['theta_grid']
        x_pcf_lo = np.hstack((0., RES['theta_grid'][:-1]))
        sumsn_100kpc = np.sum(ratio[(~np.isnan(ratio))&(x_pcf*x_conversion<0.1)])
        #if sumsn_100kpc>5:
        plt.plot(x_pcf*x_conversion, y_pcf, color=cc)

    plt.scatter(-1*c_val_norm, -1*c_val_norm, s=1, c=RES['color_grid'][sel] , cmap=pl.cm.tab20b)#, label=r"Events $1/[(R/R_c)^{\alpha}x(1+(R/R_c)^{2.2})]$")

    plt.colorbar(label=r'color VIS-y')

    plt.xlim((x_pcf.min()*x_conversion/1.5, x_pcf.max()*x_conversion*1.5))
    plt.ylim((0.1, 1000))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"$r$ [Mpc, proper] ")
    plt.ylabel(r"$\frac{N^{pairs}_{Cluster-Gal}}{N^{pairs}_{Random-Gal}}\frac{N_R}{N_D}-1$")
    plt.title(str(RES['N_D'])+' clusters, z='+str(np.round(RES['z_bar'],3)) +'\n'+r' $\Sigma [S/N(r<100kpc)] > $'+str(SN_min))
    plt.tight_layout()
    plt.savefig(p2_fig)
    plt.clf()
    print(p2_fig)




