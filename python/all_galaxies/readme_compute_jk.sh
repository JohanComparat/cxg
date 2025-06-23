#!/bin/bash
export OMP_NUM_THREADS=32
cd ~/software/cxg/python/all_galaxies

# Galaxy autocorrelation JK in the correct redshift bin

nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 4 > log_1025_jk_4.log &

nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 8 > log_1025_jk_8.log &

nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 4 > log_1075_jk_4.log &

nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 8 > log_1075_jk_8.log &

nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"BC_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"BC_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 8 > log_1025bc_jk_8.log &

nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"BC_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"BC_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 8 > log_1075bc_jk_8.log &

nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"RS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"RS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"LS10_VLIM_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02" 0.1 0.2 8 > log_1025rs_jk_8.log &

nohup python compute_wprp_Files_MaskHpx_MstarDEP_z0102_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"RS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"RS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"LS10_VLIM_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03" 0.1 0.3 8 > log_1075rs_jk_8.log &

# Cluster only


nohup python compute_wprp_Files_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_wprp-pimax100-bin0p1-01z02" \
 0.1 0.2 4 "C0xC0" 0.1 > log_C0_jk_4.log &

nohup python compute_wprp_Files_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_wprp-pimax100-bin0p1-01z02" \
 0.1 0.3 4 "C1xC1" 0.1 > log_C1_jk_4.log &



# Galaxy-Cluster cross-correlation JK in the correct redshift bin

export OMP_NUM_THREADS=32
cd ~/software/cxg/python/all_galaxies

nohup python compute_wprp_Files_MaskHpx_MstarDEP_cross_zSel_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask-01z02" \
10.25 0.1 0.2 8 "C0xG1025" 0.1 > log_C0xG1025_8.log &

nohup python compute_wprp_Files_MaskHpx_MstarDEP_cross_zSel_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"BC_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"BC_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask-01z02" \
10.25 0.1 0.2 8 "C0xG1025" 0.1 > log_C0xG1025_8_BC.log &

nohup python compute_wprp_Files_MaskHpx_MstarDEP_cross_zSel_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"RS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"RS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask-01z02" \
10.25 0.1 0.2 8 "C0xG1025" 0.1 > log_C0xG1025_8_RS.log &


export OMP_NUM_THREADS=32
cd ~/software/cxg/python/all_galaxies

python compute_wprp_Files_MaskHpx_MstarDEP_cross_zSel_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask-01z03" \
10.75 0.1 0.3 8 "C1xG1075" 0.1 > log_C1xG1075_8.log &

python compute_wprp_Files_MaskHpx_MstarDEP_cross_zSel_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"BC_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"BC_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask-01z03" \
10.75 0.1 0.3 8 "C1xG1075" 0.1 > log_C1xG1075_8_BC.log &

python compute_wprp_Files_MaskHpx_MstarDEP_cross_zSel_TrueJK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"RS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"RS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask-01z03" \
10.75 0.1 0.3 8 "C1xG1075" 0.1 > log_C1xG1075_8_RS.log &

cd ~/software/cxg/data/
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100 .
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100 .
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100 .
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100 .
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100 .
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100 .
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/C0xG1025 .
rsync -avz /data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/C1xG1075 .

git add _ANY_10.*
git add _*
git add C*
git commit -m"pcf"
git push

python plot_JK_wprp.py
