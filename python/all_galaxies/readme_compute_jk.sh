#!/bin/bash
export OMP_NUM_THREADS=32
cd ~/software/cxg/python/all_galaxies

Galaxy autocorrelation JK in the correct redshift bin

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

Galaxy-Cluster cross-correlation JK in the correct redshift bin

export OMP_NUM_THREADS=32
cd ~/software/cxg/python/all_galaxies

python compute_wprp_Files_MaskHpx_MstarDEP_cross_zSel_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask-01z02" \
10.25 0.1 0.2 4 "C0xG1025"
