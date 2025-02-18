
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, vstack
from glob import glob
import xml.etree.ElementTree as ET
import os, sys
from os.path import isfile, isdir
import pdb

#from filter_gals import filter_gals, make_unif_mask
## import from modules in VisualToolBox
#sys.path.append('modules/')
#from utils import dm_free_metadata_reader

def read_euclid_xml(xml_file, debug=False):
    """
    Simple reading of Euclid xml files
    

    Returns
    -------
    None.

    """
    tree = ET.parse(xml_file) 
    root = tree.getroot() 
    #tindex = root.find('Data/TileIndex').text
    
    tileids = []
    filenames  =  []
    for child in root:
         if debug: print(child.tag)
         if debug: print("------")
         for subchild in child:
             if debug: print("* "+ subchild.tag,subchild.text)
             if subchild.tag == "TileIndex":
                 tileids.append(subchild.text)
             if subchild.tag == "TileIndexList":
                 tileids.append(subchild.text)
             for ss in subchild:
                 if debug: print("** "+ss.tag,ss.attrib,ss.text)
                 for sss in ss:
                     if debug:  print("*** "+sss.tag,sss.attrib,sss.text)
                     if sss.tag == 'FileName':
                         filenames.append(sss.text)
                     for ssss in sss:
                         if debug:  print("*** "+ssss.tag,ssss.attrib,ssss.text)
                         if ssss.tag == 'FileName':
                             filenames.append(ssss.text)
    return tileids, filenames, tree, root

keep_fields =[
        "VIS_DET",
        "DET_QUALITY_FLAG",
        "FLUX_DETECTION_TOTAL",
        "FLUXERR_DETECTION_TOTAL",

        "OBJECT_ID",
        "RIGHT_ASCENSION",
        "DECLINATION",
        "FLUX_VIS_2FWHM_APER",
        "FLUX_Y_2FWHM_APER" ,
        "FLUX_J_2FWHM_APER" ,
        "FLUX_H_2FWHM_APER",
        "FLUX_U_EXT_DECAM_2FWHM_APER" ,
        "FLUX_G_EXT_DECAM_2FWHM_APER" ,
        "FLUX_R_EXT_DECAM_2FWHM_APER" ,
        "FLUX_I_EXT_DECAM_2FWHM_APER",
        "FLUX_Z_EXT_DECAM_2FWHM_APER",
        "FLUX_VIS_2FWHM_APER",
        "FLUXERR_Y_2FWHM_APER" ,
        "FLUXERR_J_2FWHM_APER" ,
        "FLUXERR_H_2FWHM_APER",
        "FLUXERR_U_EXT_DECAM_2FWHM_APER" ,
        "FLUXERR_G_EXT_DECAM_2FWHM_APER" ,
        "FLUXERR_R_EXT_DECAM_2FWHM_APER" ,
        "FLUXERR_I_EXT_DECAM_2FWHM_APER",
        "FLUXERR_Z_EXT_DECAM_2FWHM_APER",
        "FLAG_VIS",
        "FLAG_Y" ,
        "FLAG_J" ,
        "FLAG_H",
        "FLAG_U_EXT_DECAM" ,
        "FLAG_G_EXT_DECAM",
        "FLAG_R_EXT_DECAM" ,
        "FLAG_I_EXT_DECAM" ,
        "FLAG_Z_EXT_DECAM" ,
        ]


dir_files='/sps/euclid/OU-LE3/VMPZ-ID/newF006/sandbox/'
#dir_files='/sps/euclid/OU-LE3/VMPZ-ID/newF006/'
dir_files='/sps/euclid/OU-LE3/VMPZ-ID/REGREPROC1_R2/'
dir_files='/sps/euclid/OU-LE3/VMPZ-ID/Q1/'
phz_str = 'DpdPhzPfOutputCatalog*xml'
mer_str = 'DpdMerFinalCatalog-MER*xml'
spe_str = 'DpdSpePfOutputCatalog*xml'
subdir_results='CatRed_CL/'

mer_xmls   = glob(dir_files+mer_str)
mer_tiles = []
mer_files = []
all_t = []
for xmlf in mer_xmls:
    tid,fid, tree, root = read_euclid_xml(xmlf)
    mer_tiles.append(tid[0])
    mer_files.append(fid[0])
    t = Table.read(os.path.join(dir_files,fid[0]))
    #t = fits.open(os.path.join(dir_files,fid[0]))[1].data
    keep = (t["VIS_DET"]==1)&(t["DET_QUALITY_FLAG"]==0) & (t["FLUX_DETECTION_TOTAL"]/t["FLUXERR_DETECTION_TOTAL"]>=5) & (t["FLAG_VIS"]==0) & (t["FLAG_Y"]==0)&(t["FLAG_G_EXT_DECAM"]==0)&(t["FLAG_R_EXT_DECAM"]==0) & (t["FLUX_Y_2FWHM_APER"]/t["FLUXERR_Y_2FWHM_APER"]>=5) & (t["FLUX_G_EXT_DECAM_2FWHM_APER"]/t["FLUXERR_G_EXT_DECAM_2FWHM_APER"]>=5) & (t["FLUX_R_EXT_DECAM_2FWHM_APER"]/t["FLUXERR_R_EXT_DECAM_2FWHM_APER"]>=5) & ( -2.5*np.log10(t["FLUXERR_G_EXT_DECAM_2FWHM_APER"]/t["FLUXERR_R_EXT_DECAM_2FWHM_APER"]) > -0.5 ) & ( -2.5*np.log10(t["FLUX_VIS_2FWHM_APER"]/t["FLUX_Y_2FWHM_APER"]) > -0.5 )  & ( -2.5*np.log10(t["FLUXERR_G_EXT_DECAM_2FWHM_APER"]/t["FLUXERR_R_EXT_DECAM_2FWHM_APER"]) < 3 ) & ( -2.5*np.log10(t["FLUX_VIS_2FWHM_APER"]/t["FLUX_Y_2FWHM_APER"]) < 3 )
    t = t[keep]
    t['gr']=-2.5*np.log10(t["FLUXERR_G_EXT_DECAM_2FWHM_APER"]/t["FLUXERR_R_EXT_DECAM_2FWHM_APER"])
    t['vy']=-2.5*np.log10(t["FLUX_VIS_2FWHM_APER"]/t["FLUX_Y_2FWHM_APER"])
    t.keep_columns(["RIGHT_ASCENSION","DECLINATION",'gr','vy'])
    print(os.path.join(dir_files,fid[0]), len(t))
    if len(t)>0:
        all_t.append(t)

t_out = vstack(all_t)
t_out.write('~/color_catalog.fits', overwrite = True)
#phz_xmls = glob(dir_files+phz_str)
#phz_tiles = []
#phz_files = []
#for xmlf in phz_xmls[:1]:
    #tid,fid, tree, root = read_euclid_xml(xmlf)
    #phz_tiles.append(tid[0])
    #phz_files.append(fid[0])
    #tz = Table.read(os.path.join(dir_files,fid[0]))



#
# what is in the xml file ?
# a big mess with at least 4 nested levels :
#
for cc in root:
     print (cc)
     for c_i in cc:
         print(c_i.tag, c_i.text)
         for c_j in c_i:
             print(c_j.tag, c_j.text)
             for c_k in c_j:
                 print(c_k.tag, c_k.text)
