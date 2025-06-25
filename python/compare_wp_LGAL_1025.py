print('+'*100)
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import os
import glob
import numpy as np
from astropy.table import Table
colors = ["#67E568","#FFF000","#FFB62B","#E56124",
		  "#E53E30",
		  "#7F2353","#F911FF","#9F8CA6","#257F27","#08420D"]
syst_err = 0.02

fig_dir  ='../figures/uchuu'

os.environ['GIT_STMOD_DATA'] = os.path.join(os.environ['USERPROFILE'], "Documents\Shared\software\st_mod_data") # visible in this process + all children

os.environ['LSDR10'] = os.path.join(os.environ['USERPROFILE'], "Documents\Shared\data\legacysurvey\dr10")
LS10_bgs_vlim_dir = os.path.join(os.environ['LSDR10'], 'sweep\BGS_VLIM_Mstar')

d1025 = '../data/_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100'
d1025_BC = '../DATA/_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100'
d1025_RS = '../DATA/_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100'
dC0xG1025 = '../data/C0xG1025'
dC0xC0 = '../data/C0xC0'

lgal_dir_z0p2 = '../data/lgal'

LGAL = {}
LGAL["L0520_43.9_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.9_BC_1023_JK_100.fits" ))
LGAL["L0520_43.9_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.9_Ms_1023_JK_100.fits" ))
LGAL["L0520_43.9_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.9_RS_1023_JK_100.fits" ))

#for lmin in np.arange(42.7, 44.4,0.1):
lmin = 43.9
# for lmin in np.arange(43.9, 43.91,0.1):
l_str = str(np.round(lmin,1))
# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-RSBC-S1-clusters-LGAL.png')
plt.figure(10, (5,5))

pcf_list_1025 = np.array(glob.glob(os.path.join(dC0xG1025, '*_ANY_*NSIDE_08*.fits')))
all_wprp = []
for el in pcf_list_1025:
    all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
f_err = all_wprp.std(axis=0) / all_wprp.mean(axis=0)
f_err_tot = (f_err ** 2 + syst_err ** 2) ** 0.5
plt.plot(t_wp['rp_mid'], all_wprp.mean(axis=0) * t_wp['rp_mid'], lw=1, ls='solid', label='Galaxies x Clusters',
         color='grey')
plt.fill_between(t_wp['rp_mid'], all_wprp.mean(axis=0) * t_wp['rp_mid'] * (1 - f_err_tot),
                 all_wprp.mean(axis=0) * t_wp['rp_mid'] * (1 + f_err_tot), color='grey', alpha=0.4)

pcf_list_1025 = np.array(glob.glob(os.path.join(dC0xG1025, '*_BC_*NSIDE_08*.fits')))
all_wprp = []
for el in pcf_list_1025:
    all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
f_err = all_wprp.std(axis=0) / all_wprp.mean(axis=0)
f_err_tot = (f_err ** 2 + syst_err ** 2) ** 0.5
plt.plot(t_wp['rp_mid'], all_wprp.mean(axis=0) * t_wp['rp_mid'], lw=1, ls='solid', color='darkblue',
         label='blue-cloud')
plt.fill_between(t_wp['rp_mid'], all_wprp.mean(axis=0) * t_wp['rp_mid'] * (1 - f_err_tot),
                 all_wprp.mean(axis=0) * t_wp['rp_mid'] * (1 + f_err_tot), color='darkblue',
                 alpha=0.4)  # , label='blue-cloud')

pcf_list_1025 = np.array(glob.glob(os.path.join(dC0xG1025, '*_RS_*NSIDE_08*.fits')))
all_wprp = []
for el in pcf_list_1025:
    all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
f_err = all_wprp.std(axis=0) / all_wprp.mean(axis=0)
f_err_tot = (f_err ** 2 + syst_err ** 2) ** 0.5
plt.plot(t_wp['rp_mid'], all_wprp.mean(axis=0) * t_wp['rp_mid'], lw=1, ls='solid', color='darkred',
         label='red-sequence')
plt.fill_between(t_wp['rp_mid'], all_wprp.mean(axis=0) * t_wp['rp_mid'] * (1 - f_err_tot),
                 all_wprp.mean(axis=0) * t_wp['rp_mid'] * (1 + f_err_tot), color='darkred',
                 alpha=0.4)  # , label='red-sequence')

t_wp = LGAL["L0520_"+l_str+"_Ms_1023"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='k', label='LGAL')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']+0.05
plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='k', alpha=0.3)
t_wp = LGAL["L0520_"+l_str+"_BC_1023"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='magenta')#, label='LGAL sSFR>-11')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']+0.05
plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='magenta', alpha=0.3)
t_wp = LGAL["L0520_"+l_str+"_RS_1023"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='orange')#, label='LGAL sSFR<-11')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']+0.05
plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='orange', alpha=0.3)

y_val = np.arange(7, 3000, 1)
plt.fill_betweenx(y=y_val, x1=60*np.ones_like(y_val), x2=15*np.ones_like(y_val), color='k' , alpha=0.3)

plt.ylim((1, 1500))
plt.xlim((0.05, 60))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$")#, $\pi_{max}=100$ Mpc$/h$")
plt.legend(loc=0, fontsize=12,ncol=1, title='C0 x G1025')#, LX>'+l_str)#, title='LS10, r<19.5, 0.05<z<0.22')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)
