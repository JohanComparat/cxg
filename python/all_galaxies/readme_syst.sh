conda activate eroconda
cd ~/software/cxg/python/all_galaxies

nohup python sys_check_rhoStar.py 32 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 > log.SC_bgs_32_dr10p1 & # DONE
nohup python sys_check_rhoStar.py 64 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 > log.SC_bgs_64_dr10p1 & # DONE
nohup python sys_check_rhoStar.py 128 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 > log.SC_bgs_128_dr10p1 & # DONE


nohup python sys_check_rhoStar.py 32 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 > log_1075_syst_32.log &
nohup python sys_check_rhoStar.py 64 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 > log_1075_syst_64.log &
ohup python sys_check_rhoStar.py 128 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 > log_1075_syst_128.log &







cd ~/software/cxg/python/all_galaxies
nohup python plot_syst.py 32 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 > log.plot.32 &
nohup python plot_syst.py 64 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 > log.plot.64 &
nohup python plot_syst.py 128 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 > log.plot.128 &

nohup python plot_syst.py 32 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 > log_1075_syst_32.log &
nohup python plot_syst.py 64 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 > log_1075_syst_64.log &
nohup python plot_syst.py 128 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 > log_1075_syst_128.log &






cd ~/software/cxg/python/figures/
git add _ANY_*_syst
git commit -m"syst"
git push

