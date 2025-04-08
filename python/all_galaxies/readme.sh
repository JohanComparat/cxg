dir_clusters='/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples'
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S3.fit         
randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits
randoms-1-0-erass1sky-hod-cutselfunc20230731_S4.fits
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S4.fit
randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit
randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits
randoms-1-0-erass1sky-hod-cutselfunc20230731_S3.fits

cd /home/comparat/software/cxg/python/all_galaxies
# AUTO CORR GAL
conda activate clustering
python compute_wprp_Files_MaskHpx_MstarDEP_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"AUTOCORR_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
 10.25 0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_RAND.fits" \
"AUTOCORR_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
 10.5 0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_DATA.fits" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_RAND.fits" \
"AUTOCORR_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
 11.0 0.05 0.35

python compute_wprp_Files_MaskHpx_MstarDEP_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_DATA.fits" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_RAND.fits" \
"AUTOCORR_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
 11.5 0.05 0.35

# AUTO CORR CLU
# S0
python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.18-wprp-pimax100-bin0p1.fits" \
0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.22-wprp-pimax100-bin0p1.fits" \
0.05 0.22

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.26-wprp-pimax100-bin0p1.fits" \
0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.35-wprp-pimax100-bin0p1.fits" \
0.05 0.35

#S1
python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.18-wprp-pimax100-bin0p1.fits" \
0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.26-wprp-pimax100-bin0p1.fits" \
0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.31-wprp-pimax100-bin0p1.fits" \
0.05 0.31

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.35-wprp-pimax100-bin0p1.fits" \
0.05 0.35

#S2
python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S2_0.05_z_0.18-wprp-pimax100-bin0p1.fits" \
0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S2_0.05_z_0.26-wprp-pimax100-bin0p1.fits" \
0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S2_0.05_z_0.35-wprp-pimax100-bin0p1.fits" \
0.05 0.35

#S3
python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S3.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S3.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S3_0.05_z_0.18-wprp-pimax100-bin0p1.fits" \
0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S3.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S3.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S3_0.05_z_0.26-wprp-pimax100-bin0p1.fits" \
0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S3.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S3.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S3_0.05_z_0.35-wprp-pimax100-bin0p1.fits" \
0.05 0.35

#S4
python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S4.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S4.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S4_0.05_z_0.18-wprp-pimax100-bin0p1.fits" \
0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S4.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S4.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S4_0.05_z_0.26-wprp-pimax100-bin0p1.fits" \
0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_CLU.py 100 \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S4.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S4.fits" \
"AUTOCORR_eRASS1_VLIM_CLUSTERS_S4_0.05_z_0.35-wprp-pimax100-bin0p1.fits" \
0.05 0.35


cd /home/comparat/software/cxg/python/all_galaxies
# AUTO CORR GAL
conda activate clustering

# S0
python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.25 0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.25 0.05 0.22


python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits" \
"LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.25 0.05 0.22


python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.5 0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.75 0.05 0.31

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_DATA.fits" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.0 0.05 0.35

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_DATA.fits" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.5 0.05 0.35


cd /home/comparat/software/cxg/python/all_galaxies
# AUTO CORR GAL
conda activate clustering

# S1
python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.25 0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.5 0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.75 0.05 0.31

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits" \
"LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.75 0.05 0.31

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_DATA.fits" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.0 0.05 0.35

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_DATA.fits" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.5 0.05 0.35

cd /home/comparat/software/cxg/python/all_galaxies
conda activate clustering

# S2
python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.25 0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.5 0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_DATA.fits" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.0 0.05 0.35

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_DATA.fits" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.5 0.05 0.35


# S3
python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S3.fits" \
"eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.25 0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S3.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S3.fits" \
"eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.5 0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_DATA.fits" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S3.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S3.fits" \
"eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.0 0.05 0.35

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_DATA.fits" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S3.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S3.fits" \
"eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.5 0.05 0.35


# S4
python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits" \
"LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S4.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S4.fits" \
"eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.25 0.05 0.18

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits" \
"LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S4.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S4.fits" \
"eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
10.5 0.05 0.26

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_DATA.fits" \
"LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S4.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S4.fits" \
"eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.0 0.05 0.35

python compute_wprp_Files_MaskHpx_MstarDEP_cross_JK.py 100 \
"/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_DATA.fits" \
"LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_RAND.fits" \
"/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S4.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S4.fits" \
"eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask_JK.fits" \
11.5 0.05 0.35


cd /home/comparat/sf_Shared/data/legacysurvey/dr10/sweep/BGS_VLIM_Mstar
rsync -avz comparat@ds43:/data36s/comparat/legacysurvey/dr10/sweep/BGS_VLIM_Mstar/*wprp*.fits .
rsync -avz comparat@ds43:/home/comparat/data/eROSITA/cluster_clustering/eRASS1_CLU_VolLimSamples/*wprp*.fits .


cd $GIT_STMOD_DATA/data/validation/validation_GAS/WPRP
rsync -avz comparat@ds52:~/st_mod_data/data/validation/validation_GAS/WPRP/* .


cd wwwDir/stuff
tar -czf test2.tar.gz CF

ls *-wprp-pimax100-bin0p1-HpxMask_JK.fits
ls *-wprp-pimax100-bin0p1.fits


# continue with
# red and blue auto-corr + xcorr


LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_DATA.fits
LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_DATA.fits
LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_DATA.fits
LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_DATA.fits
LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_DATA.fits
LS10_VLIM_ANY_11.25_Mstar_12.0_0.05_z_0.35_N_0541855_DATA.fits
LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882_DATA.fits
LS10_VLIM_ANY_9.0_Mstar_12.0_0.05_z_0.08_N_0523486_DATA.fits
LS10_VLIM_ANY_9.5_Mstar_12.0_0.05_z_0.12_N_1432502_DATA.fits


corrfunc cannot be installed locally on the laptop :

conda create -n corrfunc python=3.13
conda activate corrfunc
conda install -c conda-forge gsl
conda install -c conda-forge numpy
conda install -c conda-forge scipy
conda install -c conda-forge matplotlib
conda install -c conda-forge pip

# method 1 fails
python -m pip install Corrfunc # FAILS

# method 2 FAILS
git clone https://github.com/manodeep/Corrfunc.git
cd Corrfunc
make # FAILS
make install
python -m pip install . [--user]

make tests  # run the C tests
python -m pip install pytest
python -m pytest  # run the Python tests
