#!/bin/bash
export OMP_NUM_THREADS=32
cd ~/software/cxg/python/all_galaxies

#nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
#"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 2 > log_1025_jk_2.log &
#
#nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
#"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 4 > log_1025_jk_4.log &
#
#nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
#"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 8 > log_1025_jk_8.log &
#
#nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
#"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
#"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 16 > log_1025_jk_16.log &
#
#nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
#"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
#"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
#"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
#"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 2 > log_1075_jk_2.log &
#
#nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
#"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
#"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
#"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
#"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 4 > log_1075_jk_4.log &

# TODO
nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 8 > log_1075_jk_8.log &

# TODO
nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 16 > log_1075_jk_16.log &

cd ~/software/cxg/data/
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100 .
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100 .

python plot_JK_wprp.py

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_10.25-wprp-pimax100-bin0p05-HpxMask.fits" 10.25

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_10.5-wprp-pimax100-bin0p05-HpxMask.fits" 10.5

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_10.75-wprp-pimax100-bin0p05-HpxMask.fits" 10.75

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_11.0-wprp-pimax100-bin0p05-HpxMask.fits" 11.0

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_11.25-wprp-pimax100-bin0p05-HpxMask.fits" 11.25



python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.5-wprp-pimax100-bin0p05-HpxMask.fits" 10.5

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.75-wprp-pimax100-bin0p05-HpxMask.fits" 10.75

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_11.0-wprp-pimax100-bin0p05-HpxMask.fits" 11.0

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_11.25-wprp-pimax100-bin0p05-HpxMask.fits" 11.25



python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_RAND.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-MstarMin_10.75-wprp-pimax100-bin0p05-HpxMask.fits" 10.75

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_RAND.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-MstarMin_11.0-wprp-pimax100-bin0p05-HpxMask.fits" 11.0

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_RAND.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-MstarMin_11.25-wprp-pimax100-bin0p05-HpxMask.fits" 11.25





python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_11.0-wprp-pimax100-bin0p05-HpxMask.fits" 11.0

python compute_wprp_Files_MaskHpx_MstarDEP.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_11.25-wprp-pimax100-bin0p05-HpxMask.fits" 11.25


