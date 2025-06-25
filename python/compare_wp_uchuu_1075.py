print('+'*100)
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import os, sys
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

d1075 = '../data/_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100'
d1075_BC = '../DATA/_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100'
d1075_RS = '../DATA/_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100'
dC1xG1075 = '../data/C1xG1075'
dC1xC1 = '../data/C1xC1'

uchuu_dir_z0p14 = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/WPRP', 'z0p14')
uchuu_dir_z0p19 = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/WPRP', 'z0p19')
uchuu_dir_z0p25 = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/WPRP', 'z0p25')

def get_data(file_list):
    t = Table.read(file_list[0])
    all_wp = np.array([ Table.read(el)['wprp'] for el in file_list ])
    t['wprp_JK_mean'] = np.mean(all_wp, axis=0)
    t['wprp_JK_std'] = np.std(all_wp, axis=0)
    return t

U0p14={}
U0p14["BC_10.25_LX42.8"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_BC_HaloLX_gt_42.8_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p14["RS_10.25_LX42.8"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_RS_HaloLX_gt_42.8_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p14["ANY_10.25_LX42.8"] = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_HaloLX_gt_42.8_replication_*-wprp-pimax100-2pcf-bin0p1.fits"    )))
U0p14["BC_10.25_LX42.9"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_BC_HaloLX_gt_42.9_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p14["RS_10.25_LX42.9"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_RS_HaloLX_gt_42.9_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p14["ANY_10.25_LX42.9"] = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_HaloLX_gt_42.9_replication_*-wprp-pimax100-2pcf-bin0p1.fits"    )))

U0p19={}
U0p19[ "BC_10.75_LX43.1"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_BC_HaloLX_gt_43.1_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p19[ "RS_10.75_LX43.1"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_RS_HaloLX_gt_43.1_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p19["ANY_10.75_LX43.1"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_HaloLX_gt_43.1_replication_*-wprp-pimax100-2pcf-bin0p1.fits"    )))
U0p19[ "BC_10.75_LX43.2"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_BC_HaloLX_gt_43.2_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p19[ "RS_10.75_LX43.2"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_RS_HaloLX_gt_43.2_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p19["ANY_10.75_LX43.2"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_HaloLX_gt_43.2_replication_*-wprp-pimax100-2pcf-bin0p1.fits"    )))


#
#
# S1 10.75
#
#


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS-S1-clusters.png')
plt.figure(10, (5,5))


pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_ANY_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
f_err_tot = (f_err**2 + syst_err**2)**0.5
# f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'],all_wprp.mean(axis=0)*t_wp['rp_mid'], lw=1,  ls='solid', label='Galaxies x Clusters', color='grey')
plt.fill_between(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid']*(1-f_err_tot), all_wprp.mean(axis=0)*t_wp['rp_mid']*(1+f_err_tot), color='grey', alpha=0.4)

pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_BC_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
f_err_tot = (f_err**2 + syst_err**2)**0.5
# f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid'], lw=1,  ls='solid', color='darkblue', label='blue-cloud')
plt.fill_between(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid']*(1-f_err_tot), all_wprp.mean(axis=0)*t_wp['rp_mid']*(1+f_err_tot), color='darkblue', alpha=0.4)#, label='blue-cloud')

pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_RS_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
f_err_tot = (f_err**2 + syst_err**2)**0.5
# f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid'], lw=1,  ls='solid', color='darkred', label='red-sequence')
plt.fill_between(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid']*(1-f_err_tot), all_wprp.mean(axis=0)*t_wp['rp_mid']*(1+f_err_tot), color='darkred', alpha=0.4)#, label='red-sequence')

t_wp = U0p19["ANY_10.75_LX43.2"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='k', label='UM')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='k', alpha=0.3)
t_wp = U0p19["BC_10.75_LX43.2"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='magenta')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='magenta', alpha=0.3)
t_wp = U0p19["RS_10.75_LX43.2"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='orange')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='orange', alpha=0.3)

y_val = np.arange(7, 3000, 1)
plt.fill_betweenx(y=y_val, x1=60*np.ones_like(y_val), x2=15*np.ones_like(y_val), color='k' , alpha=0.3)

plt.ylim((1, 1500))
plt.xlim((0.05, 60))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$")#, $\pi_{max}=100$ Mpc$/h$")
plt.legend(loc=0, fontsize=12,ncol=1, title='C1 x G1075')#, title='LS10, r<19.5, 0.05<z<0.22')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)
