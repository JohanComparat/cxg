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
from scipy.stats import norm
import astropy.units as u
import astropy.constants as cc
from astropy.cosmology import FlatLambdaCDM
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT

os.environ['GIT_STMOD_DATA'] = os.path.join(os.environ['USERPROFILE'], "Documents\Shared\software\st_mod_data") # visible in this process + all children

RS_model = Table.read( os.path.join( os.environ['GIT_STMOD_DATA'], 'data', 'models', 'model_GAL', 'legacy_dr10_south_v0.3_grz_z_cal_zspec_redgals_model.fit') )
z_RS = np.hstack(( RS_model['nodes'][0] ))
gr_RS = np.hstack(( np.array(RS_model['meancol']).T[0].T[0] ))
gr_RS_sigma = np.hstack(( np.array(RS_model['meancol_scatter']).T[0].T[0] ))
RS_interp = interp1d(z_RS, gr_RS)
RS_interp_up = interp1d(z_RS, gr_RS + gr_RS_sigma)
RS_interp_lo = interp1d(z_RS, gr_RS - gr_RS_sigma)
RS_interp_lsigma = interp1d(z_RS, gr_RS_sigma)

fig_dir  ='../../figures/'
dat_dir  ='../../data/'
deg_to_rad = np.pi/180.

cl_sample = sys.argv[1] # 'S0'
p2_all_dat  = np.array(glob.glob(os.path.join(dat_dir, 'Counts_DATA_'+cl_sample+'_RAND_gr_*_z_*.RSdata.npy')))

all_z         = []
all_colors    = []
all_sn_1Mpc   = []
all_sn_100kpc = []
intWTH_100kpc    = []
intWTH_100kpc_up = []
intWTH_100kpc_lo = []
intWTH_1Mpc    = []
intWTH_1Mpc_up = []
intWTH_1Mpc_lo = []
color_step = 0.05

RES = {}
RES_MAX = []
RES_MAX_100kpc = []
params_SN = []
params_WT = []
params_WT_err = []
for jj, p2_data in enumerate(p2_all_dat[:500]):
    RES[jj] = np.load(p2_data, allow_pickle='TRUE').item()
    all_z         .append(RES[jj]['z_bar'] * np.ones_like(RES[jj]['color_grid']))
    all_colors    .append(RES[jj]['color_grid']+color_step/2.)
    all_sn_1Mpc   .append(RES[jj]['sumsn_1Mpc_all'])
    all_sn_100kpc .append(RES[jj]['sumsn_100kpc_all'])
    #
    intWTH_100kpc    .append( RES[jj]['intWTH_100kpc']    )
    intWTH_100kpc_up .append( RES[jj]['intWTH_100kpc_up'] )
    intWTH_100kpc_lo .append( RES[jj]['intWTH_100kpc_lo'] )
    intWTH_1Mpc    .append( RES[jj]['intWTH_1Mpc']      )
    intWTH_1Mpc_up .append( RES[jj]['intWTH_1Mpc_up']   )
    intWTH_1Mpc_lo .append( RES[jj]['intWTH_1Mpc_lo']   )
    #
    i_max1 = np.argmax(RES[jj]['sumsn_1Mpc_all'])
    RES_MAX.append([RES[jj]['z_bar'], RES[jj]['color_grid'][i_max1]+color_step/2.])
    i_max = np.argmax(RES[jj]['sumsn_100kpc_all'])
    RES_MAX_100kpc.append([RES[jj]['z_bar'], RES[jj]['color_grid'][i_max]+color_step/2.])
    # individual figure
    p2_fig  = os.path.join(fig_dir, 'RS_gr_z_SN_'+os.path.basename(p2_data)+'.png')
    plt.figure(11, (7.4,7.5))
    #yerr = RES[jj]['sumsn_100kpc_all']*0.5*(RES[jj]['N_data_100kpc']**(-1)+RES[jj]['N_rand_100kpc']**(-1))**0.5
    #plt.errorbar(RES[jj]['color_grid']+color_step/2., RES[jj]['sumsn_100kpc_all'], xerr=color_step/2., yerr=yerr, label=r'$R_{max}=200kpc$', ls='', lw=1, c='darkred')
    yerr = RES[jj]['sumsn_100kpc_all']*0.5*(RES[jj]['N_data_100kpc']**(-1)+RES[jj]['N_rand_100kpc']**(-1))**0.5
    plt.errorbar(RES[jj]['color_grid']+color_step/2., RES[jj]['sumsn_100kpc_all'], xerr=color_step/2., yerr=yerr, label=r'$R_{max}=200kpc$', ls='', lw=1, c='darkred' )
    #plt.axvline(RS_interp(RES[jj]['z_bar']), color='k', lw=2, ls='--', label='model Kluge 24, g-r='+str(np.round(RS_interp(RES[jj]['z_bar']),3)))
    #plt.axvline(RS_interp_up(RES[jj]['z_bar']), color='k', lw=2, ls=':')
    #plt.axvline(RS_interp_lo(RES[jj]['z_bar']), color='k', lw=2, ls=':')
    #plt.axvline(RES[jj]['color_grid'][i_max]+color_step/2., label=r'$R_{max}=200kpc$, g-r=' +str(np.round(RES[jj]['color_grid'][i_max]+color_step/2.,3)), c='darkred' , lw=2, ls='--')
    #plt.axvline(RES[jj]['color_grid'][i_max1]+color_step/2., label=r'$R_{max}=500kpc$, g-r='+str(np.round(RES[jj]['color_grid'][i_max1]+color_step/2.,3)), c='darkgreen' , lw=2, ls='--')
    x_values = np.arange(0,2,0.001) #RES[jj]['color_grid'].min()+color_step/2.
    x_loc = RES[jj]['color_grid'][i_max1]+color_step/2.
    scale = RS_interp_lsigma(RES[jj]['z_bar'])
    y_values = norm.pdf(x_values, loc=RS_interp(RES[jj]['z_bar']), scale=scale)
    plt.plot(x_values, y_values*np.max(RES[jj]['sumsn_100kpc_all'])/y_values.max(), label='model Kluge 24, g-r='+str(np.round(RS_interp(RES[jj]['z_bar']),3))+r', $\sigma=$'+str(np.round(scale,3)), color='k')
    #scale = 0.1
    #y_values = norm.pdf(x_values, loc=x_loc, scale=scale)
    #plt.plot(x_values, y_values*np.max(RES[jj]['sumsn_100kpc_all'])/y_values.max(), label=r'$\sigma=$'+str(np.round(scale,2)))
    #scale = 0.2
    #y_values = norm.pdf(x_values, loc=x_loc, scale=scale)
    #plt.plot(x_values, y_values*np.max(RES[jj]['sumsn_100kpc_all'])/y_values.max(), label=r'$\sigma=$'+str(np.round(scale,2)))
    fun = lambda c_val, c_0, c_m, c_scale, c_loc : ( norm.pdf(c_val, loc=c_loc, scale=c_scale) + c_0 ) * c_m
    p_out, c_out =curve_fit(fun, RES[jj]['color_grid']+color_step/2., RES[jj]['sumsn_100kpc_all'], sigma = 10 ,p0=( 0, np.max(RES[jj]['sumsn_100kpc_all']), RS_interp_lsigma(RES[jj]['z_bar']), RS_interp(RES[jj]['z_bar']) ))
    # + (RES[jj]['intWTH_100kpc_up']-RES[jj]['intWTH_100kpc_lo'])/10
    y_values = fun(x_values, p_out[0], p_out[1], p_out[2], p_out[3])
    plt.plot(x_values, y_values*np.max(RES[jj]['sumsn_100kpc_all'])/y_values.max(), label=r'fit g-r='+str(np.round(p_out[3],3)) + ', $\sigma=$'+str(np.round(p_out[2],3)), color='r' )
    plt.xlabel(r'color $g-r$')
    plt.ylabel(r'$\Sigma [S/N(r<R_{max})]$')
    #plt.yscale('log')
    plt.xlim((0,2))
    plt.ylim((1,2*np.max(RES[jj]['sumsn_100kpc_all'])))
    plt.title('clusters '+cl_sample +',z='+str(np.round(RES[jj]['z_bar'],2)))
    plt.legend(loc=0, fontsize=10, framealpha=0.6)
    plt.tight_layout()
    plt.savefig(p2_fig)
    plt.clf()
    print(p2_fig)
    params_SN.append(np.hstack((p_out, RES[jj]['z_bar'])))

    p2_fig  = os.path.join(fig_dir, 'RS_gr_z_intWTH_'+os.path.basename(p2_data)+'-merged.png')
    #plt.figure(15, (5.5, 5.5))
    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        sharex=True,
        gridspec_kw={'height_ratios': [4, 1]},
        figsize=(5.5, 9)
    )
    s1 = (RES[jj]['intWTH_1Mpc']!=np.inf)&(RES[jj]['intWTH_1Mpc']>0.0)&(np.isinf(RES[jj]['intWTH_1Mpc_up'])==False) & (np.isinf(RES[jj]['intWTH_1Mpc_lo'])==False)
    i_max2 = np.argmax(RES[jj]['intWTH_1Mpc'][s1])
    x_loc = RES[jj]['color_grid'][s1][i_max2]+color_step/2.
    #i_max3 = np.argmax(RES[jj]['intWTH_100kpc'])
    #ax1.plot(RES[jj]['color_grid']+color_step/2., RES[jj]['intWTH_100kpc'], ls='', color='darkgreen', marker='+')#, xerr=color_step/2., yerr=[RES[jj]['intWTH_100kpc']-RES[jj]['intWTH_100kpc_lo'], RES[jj]['intWTH_100kpc_up']-RES[jj]['intWTH_100kpc']], label=r'$R_{max}=200kpc$', ls='', lw=1, c='darkgreen')
    #ax1.fill_between(RES[jj]['color_grid']+color_step/2., RES[jj]['intWTH_100kpc_lo'], RES[jj]['intWTH_100kpc_up'], color='darkgreen', alpha=0.3)
    ax1.plot(RES[jj]['color_grid'][s1]+color_step/2., RES[jj]['intWTH_1Mpc'][s1], label=r'DATA ', ls='', marker='x', color='darkred') #  , xerr=color_step/2., yerr=[RES[jj]['intWTH_1Mpc'][s1]-RES[jj]['intWTH_1Mpc_lo'][s1], RES[jj]['intWTH_1Mpc_up'][s1]-RES[jj]['intWTH_1Mpc'][s1]])
    ax1.fill_between(RES[jj]['color_grid'][s1]+color_step/2., RES[jj]['intWTH_1Mpc_lo'][s1], RES[jj]['intWTH_1Mpc_up'][s1], color='darkred', alpha=0.3)
    #ax1.axvline(x_loc, label=r'$R_{max}=500kpc$, g-r=' +str(np.round(x_loc,3)), c='darkred' , ls='--', lw=2)
    #ax1.axvline(RES[jj]['color_grid'][i_max3]+color_step/2., label=r'$R_{max}=200kpc$, g-r='+str(np.round(RES[jj]['color_grid'][i_max3]+color_step/2.,3)), c='darkgreen' , ls='--', lw=2)
    x_values = np.arange(0,2,0.001) #RES[jj]['color_grid'].min()+color_step/2.
    scale = RS_interp_lsigma(RES[jj]['z_bar'])
    y_values = norm.pdf(x_values, loc=RS_interp(RES[jj]['z_bar']), scale=scale)
    ax1.plot(x_values, y_values*np.max(RES[jj]['intWTH_1Mpc'][s1])/y_values.max(), color='k', label='model Kluge 24, g-r='+str(np.round(RS_interp(RES[jj]['z_bar']),3)) + r', $\sigma=$'+str(np.round(scale,3)))
    #ax1.axvline(RS_interp(RES[jj]['z_bar']), color='k', lw=2, ls='--',  )
    #ax1.axvline(RS_interp_up(RES[jj]['z_bar']), color='k', lw=2, ls=':')
    #ax1.axvline(RS_interp_lo(RES[jj]['z_bar']), color='k', lw=2, ls=':')
    fun = lambda c_val, c_0, c_m, c_scale, c_loc : ( norm.pdf(c_val, loc=c_loc, scale=c_scale) + c_0 ) * c_m
    p_out, c_out =curve_fit(fun, RES[jj]['color_grid'][s1]+color_step/2., RES[jj]['intWTH_1Mpc'][s1] #sigma = (RES[jj]['intWTH_1Mpc_up'][s1]-RES[jj]['intWTH_1Mpc_lo'][s1])/10,.
                            ,p0=( 0, np.max(RES[jj]['intWTH_1Mpc'][s1]), RS_interp_lsigma(RES[jj]['z_bar']), RS_interp(RES[jj]['z_bar']) ))
    y_values = fun(x_values, p_out[0], p_out[1], p_out[2], p_out[3])
    ax1.plot(x_values, y_values*np.max(RES[jj]['intWTH_1Mpc'][s1])/y_values.max(), label=r'fit g-r='+str(np.round(p_out[3],3)) + ', $\sigma=$'+str(np.round(p_out[2],3)), color='r' )
    #scale = 0.2
    #y_values = norm.pdf(x_values, loc=x_loc, scale=scale)
    #ax1.plot(x_values, y_values*np.max(RES[jj]['intWTH_1Mpc'][s1])/y_values.max(), label=r'$\sigma=$'+str(np.round(scale,2)))
    #ax1.set_xlabel(r'color $g-r$')
    ax1.set_xlim((RS_interp(RES[jj]['z_bar'])-0.5,RS_interp(RES[jj]['z_bar'])+0.5))
    ax1.set_ylabel(r'$w_{500kpc}$')
    #ax1.set_yscale('log')
    ax1.set_title('clusters '+cl_sample +',z='+str(np.round(RES[jj]['z_bar'],3)))
    ax1.legend(loc=0, fontsize=10, framealpha=0.6)
    #ax1.tight_layout()
    #plt.savefig(p2_fig)
    #plt.clf()
    print(p2_fig)
    params_WT.append(np.hstack((p_out, RES[jj]['z_bar'])))
    params_WT_err.append([ c_out[3][3]**0.5, c_out[2][2]**0.5 ])
    #p2_fig  = os.path.join(fig_dir, 'RS_gr_z_intWTH_'+os.path.basename(p2_data)+'-residual.png')
    #plt.figure(8, (5.6, 3.5))
    ax2.fill_between(RES[jj]['color_grid'][s1]+color_step/2., RES[jj]['intWTH_1Mpc_lo'][s1]/RES[jj]['intWTH_1Mpc'][s1], RES[jj]['intWTH_1Mpc_up'][s1]/RES[jj]['intWTH_1Mpc'][s1], color='darkred', alpha=0.3, label='Relative uncertainties')
    y_model = fun(RES[jj]['color_grid'][s1]+color_step/2., p_out[0], p_out[1], p_out[2], p_out[3])
    ax2.plot(RES[jj]['color_grid'][s1]+color_step/2., RES[jj]['intWTH_1Mpc'][s1]/y_model, label=r'DATA / model', ls='--', marker='o', color='darkred')
    y_val_MK24 = norm.pdf(RES[jj]['color_grid'][s1]+color_step/2., loc=RS_interp(RES[jj]['z_bar']), scale=RS_interp_lsigma(RES[jj]['z_bar']))
    y_MK24 = y_val_MK24*np.max(RES[jj]['intWTH_1Mpc'][s1])/y_val_MK24.max()
    ax2.plot(RES[jj]['color_grid'][s1]+color_step/2., RES[jj]['intWTH_1Mpc'][s1]/y_MK24, label=r'DATA / model K24', ls='--', marker='s', color='k')
    ax2.axvline(p_out[3], color='darkred', lw=2, ls='--',  )
    ax2.axvline(p_out[3]+p_out[2], color='darkred', lw=2, ls=':')
    ax2.axvline(p_out[3]-p_out[2], color='darkred', lw=2, ls=':')
    ax2.set_xlabel(r'color $g-r$')
    ax2.set_xlim((RS_interp(RES[jj]['z_bar'])-0.5,RS_interp(RES[jj]['z_bar'])+0.5))
    ax2.set_ylabel(r'data/model')
    ax2.set_ylim((0,2))
    #ax2.set_yscale('log')
    #ax2.set_title('clusters '+cl_sample +',z='+str(np.round(RES[jj]['z_bar'],3)))
    #ax2.legend(loc=3, fontsize=10, framealpha=0.6)
    plt.tight_layout()
    plt.savefig(p2_fig)
    plt.clf()
    print(p2_fig)


params_WT = np.transpose(params_WT)
params_WT_err = np.transpose(params_WT_err)
params_SN = np.transpose(params_SN)

RES_MAX = np.transpose(RES_MAX)
RES_MAX_100kpc = np.transpose(RES_MAX_100kpc)

all_z         = np.hstack((all_z        ))
all_colors    = np.hstack((all_colors   ))
all_sn_1Mpc   = np.hstack((all_sn_1Mpc  ))
all_sn_100kpc = np.hstack((all_sn_100kpc))
intWTH_100kpc     = np.hstack((  intWTH_100kpc    ))
intWTH_100kpc_up  = np.hstack((  intWTH_100kpc_up ))
intWTH_100kpc_lo  = np.hstack((  intWTH_100kpc_lo ))
intWTH_1Mpc     = np.hstack((  intWTH_1Mpc    ))
intWTH_1Mpc_up  = np.hstack((  intWTH_1Mpc_up ))
intWTH_1Mpc_lo  = np.hstack((  intWTH_1Mpc_lo ))





all_sn_1Mpc[all_sn_1Mpc<=0]=0.001
all_sn_100kpc[all_sn_100kpc<=0]=0.001
RS = {}
RS['all_sn_100kpc'] = all_sn_100kpc
RS['all_sn_1Mpc']   = all_sn_1Mpc
RS['all_colors']    = all_colors
RS['all_z']         = all_z
RS['intWTH_100kpc']    = intWTH_100kpc
RS['intWTH_100kpc_up'] = intWTH_100kpc_up
RS['intWTH_100kpc_lo'] = intWTH_100kpc_lo
RS['intWTH_1Mpc']    = intWTH_1Mpc
RS['intWTH_1Mpc_up'] = intWTH_1Mpc_up
RS['intWTH_1Mpc_lo'] = intWTH_1Mpc_lo

# recast in a matrix
#mat = np.zeros(( len(np.unique(RS['all_z'])), len(np.unique(RS['all_colors'])) ))
#for ii, zz in enumerate(np.unique(RS['all_z'])):
    #for jj, yy in enumerate(np.unique(RS['all_colors'])):
        #if len(RS['intWTH_1Mpc'][(RS['all_z']==zz)&(RS['all_colors']==yy)])==1:
            #mat[ii, jj] = RS['intWTH_1Mpc'][(RS['all_z']==zz)&(RS['all_colors']==yy)]
#mat = np.transpose(mat)
# plot CT

#p2_fig  = os.path.join(fig_dir, 'RS_gr_z_WT200kpc_'+cl_sample+'.png')
#plt.figure(3, (5,8))
#ok = (RS['intWTH_100kpc']>0) & (RS['intWTH_100kpc']<np.inf)
#plt.scatter(RS['all_z'][ok], RS['all_colors'][ok], s=8, marker='s', c=RS['intWTH_100kpc'][ok] ,  label='GalxClusters', alpha=0.8)#, vmax=15) cmap=pl.cm.cool_r,
##plt.imshow(mat)#RS['all_z'], RS['all_colors'], s=10, marker='s', c=RS['intWTH_100kpc'] , cmap=pl.cm.cool_r, label='GalxClusters', alpha=0.2)#, vmax=15)
#plt.plot(z_RS, gr_RS, 'k--',lw=2,label='model Kluge 24')
#plt.plot(z_RS, gr_RS-gr_RS_sigma, 'k:', lw=2)
#plt.plot(z_RS, gr_RS+gr_RS_sigma, 'k:', lw=2)
#plt.plot(params_WT[-1], params_WT[3], 'r--', lw=2,label=r'fit, this work ',c='darkred')
#plt.plot(params_WT[-1], params_WT[3]+params_WT[2], 'r:', lw=2,c='darkred')
#plt.plot(params_WT[-1], params_WT[3]-params_WT[2], 'r:', lw=2,c='darkred')
####plt.plot(RES_MAX_100kpc[0], RES_MAX_100kpc[1], ls='--',lw=2,label='highest S/N 100kpc',c='darkgreen')
#plt.colorbar(label=r'$\int [w(r<200kpc))$')
#plt.xlim((RS['all_z'].min()-0.01, RS['all_z'].max()+0.01))
#plt.ylim((0.7, 1.9))#RS['all_colors'].min()-0.01, RS['all_colors'].max()+0.01))
#plt.ylabel(r'color $g-r$')
#plt.xlabel(r"$z$")
#plt.title('clusters '+cl_sample )
#plt.legend(loc=4, fontsize=12, framealpha=0.95)
#plt.tight_layout()
#plt.savefig(p2_fig)
#plt.clf()
#print(p2_fig)

p2_fig  = os.path.join(fig_dir, 'RS_gr_z_WT200kpc_'+cl_sample+'.png')
plt.figure(5, (5,9))
ok = (RS['intWTH_100kpc']>0) & (RS['intWTH_100kpc']<np.inf)
plt.scatter(RS['all_z'][ok], RS['all_colors'][ok], s=95, marker='s', c=np.log10(RS['intWTH_100kpc'][ok]) ,  label='C2 x BGS', alpha=0.5, linewidths=0)#, vmax=15) cmap=pl.cm.cool_r,
plt.plot(z_RS, gr_RS, 'k--',lw=2,label='model Kluge 24')
plt.plot(z_RS, gr_RS-gr_RS_sigma, 'k:', lw=2)
plt.plot(z_RS, gr_RS+gr_RS_sigma, 'k:', lw=2)
plt.plot(params_WT[-1], params_WT[3], 'r--', lw=2,label=r'fit, this work',c='darkred')
plt.plot(params_WT[-1], params_WT[3]+params_WT[2], 'r:', lw=2,c='darkred')
plt.plot(params_WT[-1], params_WT[3]-params_WT[2], 'r:', lw=2,c='darkred')
#plt.plot(RES_MAX_100kpc[0], RES_MAX_100kpc[1], ls='--',lw=2,label='highest S/N 100kpc',c='darkgreen')
plt.colorbar(label=r'$\log_{10}(w_{200kpc})$')
plt.xlim((0.15, 0.35))#RS['all_z'].min()-0.01, RS['all_z'].max()+0.01))
#plt.ylim((RS['all_colors'].min()-0.01, RS['all_colors'].max()+0.01))
plt.ylim((0.7, 1.9))#RS['all_colors'].min()-0.01, RS['all_colors'].max()+0.01))
plt.ylabel(r'color $g-r$')
plt.xlabel(r"$z$")
#plt.title('clusters '+cl_sample )
plt.legend(loc=4, fontsize=12, framealpha=0.95)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

p2_fig  = os.path.join(fig_dir, 'RS_gr_z_WT500kpc_'+cl_sample+'.png')
plt.figure(5, (5,9))
ok = (RS['intWTH_1Mpc']>0) & (RS['intWTH_1Mpc']<np.inf)
plt.scatter(RS['all_z'][ok], RS['all_colors'][ok], s=95, marker='s', c=np.log10(RS['intWTH_1Mpc'][ok]) ,  label='C2 x BGS', alpha=0.5, linewidths=0)#, vmax=15) cmap=pl.cm.cool_r,
plt.plot(z_RS, gr_RS, 'k--',lw=2,label='model Kluge 24')
plt.plot(z_RS, gr_RS-gr_RS_sigma, 'k:', lw=2)
plt.plot(z_RS, gr_RS+gr_RS_sigma, 'k:', lw=2)
plt.plot(params_WT[-1], params_WT[3], 'r--', lw=2,label=r'fit, this work',c='darkred')
plt.plot(params_WT[-1], params_WT[3]+params_WT[2], 'r:', lw=2,c='darkred')
plt.plot(params_WT[-1], params_WT[3]-params_WT[2], 'r:', lw=2,c='darkred')
#plt.plot(RES_MAX_100kpc[0], RES_MAX_100kpc[1], ls='--',lw=2,label='highest S/N 100kpc',c='darkgreen')
plt.colorbar(label=r'$\log_{10}(w_{500kpc})$')
plt.xlim((0.15, 0.35))#RS['all_z'].min()-0.01, RS['all_z'].max()+0.01))
#plt.ylim((RS['all_colors'].min()-0.01, RS['all_colors'].max()+0.01))
plt.ylim((0.7, 1.9))#RS['all_colors'].min()-0.01, RS['all_colors'].max()+0.01))
plt.ylabel(r'color $g-r$')
plt.xlabel(r"$z$")
#plt.title('clusters '+cl_sample )
plt.legend(loc=4, fontsize=12, framealpha=0.95)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

p2_fig  = os.path.join(fig_dir, 'RS_gr_z_SN200kpc_'+cl_sample+'.png')
plt.figure(11, (7.4,7.5))
plt.scatter(RS['all_z'], RS['all_colors'], s=100, marker='s', c=np.log10(RS['all_sn_100kpc']) , cmap=pl.cm.cool_r, vmin=0, label='GalxClusters', alpha=0.3)#, vmax=15)
plt.plot(z_RS, gr_RS, 'k--', lw=2, label='model Kluge 24')
plt.plot(z_RS, gr_RS-gr_RS_sigma, 'k:', lw=2)
plt.plot(z_RS, gr_RS+gr_RS_sigma, 'k:', lw=2)
plt.plot(RES_MAX[0], RES_MAX[1], ls='--',lw=2,label='highest S/N 500kpc',c='darkred')
plt.plot(RES_MAX_100kpc[0], RES_MAX_100kpc[1], ls='--',lw=2,label='highest S/N 200kpc',c='darkgreen')
plt.plot(params_WT[-1], params_WT[3], 'r--', label=r'$\int [w(r<200kpc))$')
plt.plot(params_WT[-1], params_WT[3]+params_WT[2], 'r:')
plt.plot(params_WT[-1], params_WT[3]-params_WT[2], 'r:')
plt.colorbar(label=r'log10 $\Sigma [S/N(r<200kpc)]$')
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

p2_fig  = os.path.join(fig_dir, 'RS_gr_z_SN500kpc_'+cl_sample+'.png')
plt.figure(11, (7.4,7.5))
plt.scatter(RS['all_z'], RS['all_colors'], s=100, marker='s', c=np.log10(RS['all_sn_1Mpc']) , cmap=pl.cm.cool_r, vmin=0, label='GalxClusters', alpha=0.3)#, vmax=15)
plt.plot(z_RS, gr_RS, 'k--',lw=2,label='model Kluge 24')
plt.plot(z_RS, gr_RS-gr_RS_sigma, 'k:', lw=2)
plt.plot(z_RS, gr_RS+gr_RS_sigma, 'k:', lw=2)
plt.plot(RES_MAX[0], RES_MAX[1], ls='--',lw=2,label='highest S/N 500kpc',c='darkred')
plt.plot(RES_MAX_100kpc[0], RES_MAX_100kpc[1], ls='--',lw=2,label='highest S/N 200kpc',c='darkgreen')
plt.colorbar(label=r'log10 $\Sigma [S/N(r<500kpc)]$')
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


for e0, e1, e1e, e2, e2e in zip(params_WT[-1], params_WT[3], params_WT_err[0], params_WT[2], params_WT_err[1]):
    print(np.round(e0,3), '&', np.round(e1,4), r'$\pm$', np.round(e1e,4), '&', np.round(e2,4), r'$\pm$', np.round(e2e,4), '\\\\')

