
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, vstack
from glob import glob
import xml.etree.ElementTree as ET
import os, sys
from os.path import isfile, isdir
import pdb

from filter_gals import filter_gals, make_unif_mask
# import from modules in VisualToolBox
sys.path.append('modules/')
from utils import dm_free_metadata_reader

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
    return tileids, filenames

    
def read_vmpz_xml(xml_file, debug=False):
    """
    Simple reading of Euclid xml files
    

    Returns
    -------
    None.

    """
    tree = ET.parse(xml_file) 
    root = tree.getroot()     
    
    tileids = []
    filenames  =  []
    mfilter = []
    ftype = []
    for child in root:
         if debug: print(child.tag)
         if debug: print("------")
         for subchild in child:
             if debug: print("* "+ subchild.tag,subchild.text)
             if child.tag=='Data': ftype.append(subchild.tag)
             for ss in subchild:
                 if debug: print("** "+ss.tag,ss.attrib,ss.text)
                 if ss.tag == "Filter":
                     mfilter.append(ss.text)
                 for sss in ss:
                     if debug:  print("*** "+sss.tag,sss.attrib,sss.text)
                     if sss.tag == 'FileName':
                         filenames.append(sss.text)
                     if sss.tag == 'TileIndexList':
                         tileids.append(sss.text)                        
    return ftype, mfilter, tileids, filenames



def main_merge(dir_files, phz_str, mer_str,spe_str, mask_str, add_spe=False,subdir_results="./"):
    
    """
    
    """
    print("starting merging")
    cluster_fields =[
                    "OBJECT_ID",
                    "RIGHT_ASCENSION",
                    "DECLINATION",
                    "PHZ_PDF",
                    "PHZ_MEDIAN",
                    "PHZ_MODE_1",
                    "PHZ_FLAGS",
                    "PHZ_CORRECTION",
                    "SPE_Z" ,
                    "SPE_Z_ERR" ,
                    "FLUX_VIS_APER",
                    "FLUX_Y_APER" ,
                    "FLUX_J_APER" ,
                    "FLUX_H_APER",
                    "FLUX_NIR_STACK_APER" ,
                    "FLUX_U_EXT_DECAM_APER" ,
                    "FLUX_G_EXT_DECAM_APER" ,
                    "FLUX_R_EXT_DECAM_APER" ,
                    "FLUX_I_EXT_DECAM_APER",
                    "FLUX_Z_EXT_DECAM_APER",
                    "FLUX_U_EXT_LSST_APER",
                    "FLUX_G_EXT_LSST_APER",
                    "FLUX_R_EXT_LSST_APER" ,
                    "FLUX_I_EXT_LSST_APER",
                    "FLUX_Z_EXT_LSST_APER" ,
                    "FLUX_U_EXT_MEGACAM_APER",
                    "FLUX_R_EXT_MEGACAM_APER",
                    "FLUX_G_EXT_JPCAM_APER",
                    "FLUX_I_EXT_PANSTARRS_APER",
                    "FLUX_Z_EXT_PANSTARRS_APER" ,
                    "FLUX_G_EXT_HSC_APER" ,
                    "FLUX_Z_EXT_HSC_APER",
                    "FLUX_Y_TEMPLFIT" ,
                    "FLUX_J_TEMPLFIT" ,
                    "FLUX_H_TEMPLFIT",
                    "FLUX_U_EXT_DECAM_TEMPLFIT" ,
                    "FLUX_G_EXT_DECAM_TEMPLFIT",
                    "FLUX_R_EXT_DECAM_TEMPLFIT" ,
                    "FLUX_I_EXT_DECAM_TEMPLFIT",
                    "FLUX_Z_EXT_DECAM_TEMPLFIT",
                    "FLUX_U_EXT_LSST_TEMPLFIT" ,
                    "FLUX_G_EXT_LSST_TEMPLFIT" ,
                    "FLUX_R_EXT_LSST_TEMPLFIT" ,
                    "FLUX_I_EXT_LSST_TEMPLFIT",
                    "FLUX_Z_EXT_LSST_TEMPLFIT",
                    "FLUX_U_EXT_MEGACAM_TEMPLFIT",
                    "FLUX_R_EXT_MEGACAM_TEMPLFIT",
                    "FLUX_G_EXT_JPCAM_TEMPLFIT", 
                    "FLUX_I_EXT_PANSTARRS_TEMPLFIT",
                    "FLUX_Z_EXT_PANSTARRS_TEMPLFIT" ,
                    "FLUX_G_EXT_HSC_TEMPLFIT" ,
                    "FLUX_Z_EXT_HSC_TEMPLFIT" ,
                    "FLUX_H_TOTAL_UNIF" ,
                    "FLAG_VIS",
                    "FLAG_Y" ,
                    "FLAG_J" ,
                    "FLAG_H",
                    "FLAG_NIR_STACK" ,
                    "FLAG_U_EXT_DECAM" ,
                    "FLAG_G_EXT_DECAM",
                    "FLAG_R_EXT_DECAM" ,
                    "FLAG_I_EXT_DECAM" ,
                    "FLAG_Z_EXT_DECAM" ,
                    "FLAG_U_EXT_LSST" ,
                    "FLAG_G_EXT_LSST",
                    "FLAG_R_EXT_LSST" ,
                    "FLAG_I_EXT_LSST" ,
                    "FLAG_Z_EXT_LSST",
                    "FLAG_U_EXT_MEGACAM" ,
                    "FLAG_R_EXT_MEGACAM" ,
                    "FLAG_G_EXT_JPCAM" ,
                    "FLAG_I_EXT_PANSTARRS" ,
                    "FLAG_Z_EXT_PANSTARRS" ,
                    "FLAG_G_EXT_HSC" ,
                    "FLAG_Z_EXT_HSC" ,
                    "DEBLENDED_FLAG",
                    "DET_QUALITY_FLAG",
                    "POINT_LIKE_PROB",
                    "GAL_EBV" ]

    # Add to DM10:
    cluster_fields += ["VIS_DET","POINT_LIKE_PROB", "SPURIOUS_FLAG", "SPURIOUS_PROB"]

    phz_xmls = glob(dir_files+phz_str)
    mer_xmls   = glob(dir_files+mer_str)
    spe_xmls   = glob(dir_files+spe_str)
    mask_xmls  = glob(dir_files+mask_str)
    
    #print(spe_xmls)
    phz_tiles = []
    phz_files = [] 
    for xmlf in phz_xmls: 
        tid,fid = read_euclid_xml(xmlf)
        phz_tiles.append(tid[0])
        phz_files.append(fid[0])


    mer_tiles = []
    mer_files = []
    for xmlf in mer_xmls: 
        tid,fid = read_euclid_xml(xmlf)
        mer_tiles.append(tid[0])
        mer_files.append(fid[0])
        
    spe_tiles = []
    spe_files = []
    line_files = []
    for xmlf in spe_xmls: 
        tid,fid = read_euclid_xml(xmlf)
        #print(tid,fid)
        spe_tiles.append(tid[0])
        spe_files.append(fid[0])
        line_files.append(fid[1])




    mer_tiles = np.asarray(mer_tiles)
    phz_tiles = np.asarray(phz_tiles)
    spe_tiles = np.asarray(spe_tiles)
    
    
    
    okt,mer_idx,phz_idx  =np.intersect1d(mer_tiles,phz_tiles,return_indices=True)
    
    okn,nmer_idx,spe_idx  =np.intersect1d(mer_tiles[mer_idx],spe_tiles,return_indices=True)
    okt = okt[nmer_idx]
    mer_idx = mer_idx[nmer_idx]
    phz_idx = phz_idx[nmer_idx]
    
    ## List of masks (for diffferent selected 
    tid_mask={}
    msk_tiles= []
    msk_files = []
    
    filt_bands=["NIR_Y","NIR_H","NIR_J", "VIS"]
    for xmlf in mask_xmls: 
        #tid,fid = read_euclid_xml(xmlf)
        #fid_indx=[( msk_fits_elem in f) for f in fid]
        
        filt = dm_free_metadata_reader(xmlf, "Data.CoverageMaskHealpixParams.FilterList.Filter.Name")
        tid = dm_free_metadata_reader(xmlf, "Data.CoverageMaskHealpixParams.PatchTileList.TileIndexList")
        fid = dm_free_metadata_reader(xmlf, "Data.CoverageMaskHealpix.DataContainer.FileName")
        if filt not in filt_bands:
            continue
        if tid not in tid_mask:
            tid_mask[tid] = [fid]
        else:
            tid_mask[tid].append(fid)
        
    all_mask_tids = list(tid_mask.keys())
    #print(all_mask_tids)
    #pdb.set_trace()

    """
    
    """

    print("***************************")
    indx_msk=0
    
    if not isdir(dir_files+subdir_results): os.mkdir(dir_files+subdir_results)
    stat_file = dir_files+subdir_results+'stats.txt'
    
    with open(stat_file, 'w') as the_file:
        the_file.write(f"# TileNb \t SelectedSrc \t InitialSrc  EfectiveArea[deg] \t  WeightedArea[deg] \t RelativeEffArea \n")


    for midx,pidx,sidx in zip(mer_idx[2024:],phz_idx[2024:],spe_idx[2024:]):
        
        fm = mer_files[midx]
        fp = phz_files[pidx]
        if add_spe:
            fs = spe_files[sidx]
            fl = line_files[sidx]

        ctile = str(mer_tiles[midx])
        print("*************************************")
        print("*** Processing tile ",ctile)
        print("*** Processing MER and PHZ files : ",fm, fp)
        print("************************************")

        #tile_nb=str(okt[indx_msk])
        mask_files = None
        if ctile in all_mask_tids:
            mask_files = tid_mask[ctile]
            #print("   ** Mask file : ",tid_mask[ctile])
            #if tile_nb not in fm:
            #    pdb.set_trace()
            mask_files = [dir_files + mf for mf in mask_files]
            # Keep only files that acutally IOErrorist
            mask_files= [f for f in mask_files if isfile(f)]
            #print(f"*** For tile {ctile} processing Mask files: ", mask_files)
            
        else:
            print(f" No mask tile {ctile} found")
        #indx_msk+=1
        #fm = dir_files+'data/'+fm
        #fp = dir_files+'data/'+fp
        #fs = dir_files+'data/'+fs
        
        fm = dir_files+fm
        fp = dir_files+fp
        if add_spe:
            fs = dir_files+fs
            fl = dir_files+fl
        
        #print(fm)
        #pdb.set_trace()
        if (os.path.isfile(fm) & os.path.isfile(fp)): #& os.path.isfile(fs)):
            #pdb.set_trace()
            try:
                tab_m = Table.read(fm,format='fits')
            except IOErr:
                print("Wrong MER file")
                continue
            try:
                tab_p = Table.read(fp,hdu=1,format='fits')
            except IOError:
                print("Wrong PHZ file")
                continue
            
            allids, mgal,pgal = np.intersect1d(tab_m['OBJECT_ID'].data, tab_p['OBJECT_ID'].data,return_indices=True)
            if add_spe:
                try:
                    tab_s1 = Table.read(fs,hdu=1,format='fits')
                    tab_s2 = Table.read(fs,hdu=2,format='fits')
                    tab_s3 = Table.read(fs,hdu=1,format='fits')
                    tab_l = Table.read(fl,hdu=3,format='fits')
                except IOError:
                    print("Wrong SPE file")
                    continue
                #tab_s = tab_s[tab_s['SPE_RANK'] ==0]
                allidsspe, mgaln,sgal = np.intersect1d(tab_m['OBJECT_ID'][mgal].data, tab_s['OBJECT_ID'].data,return_indices=True)
                #mgal = mgal[mgaln]
                #pgal = pgal[mgaln]
            
            #print(" *** Merging the two following lines: ")
            
            if True:
                print("OK True")
                try:
                    tab_ok = tab_m[mgal]
                    
                except IOError:
                    print("************* Unable to find mgal at pos")
                    continue

                names =['PHZ_PDF','PHZ_MEDIAN','PHZ_70_INT','PHZ_90_INT','PHZ_95_INT','PHZ_MODE_1','PHZ_MODE_1_AREA','PHZ_FLAGS','BIAS_ID','PHOTOMETRIC_SYSTEM','PHZ_CLASSIFICATION','PHZ_WEIGHT']
                # + for shear: 
                # names +=['PHZ_MODE_2','PHZ_MODE_2_AREA,'TOM_BIN_ID','ALT_TOM_BIN_ID', 'PHZ_WEIGHT']
                for name in names:
                    tab_ok.add_column(tab_p[pgal][name],copy=True)

                # ADDED: truncated PHZ:
                tab_ok["PHZ_PDF"]=tab_ok["PHZ_PDF"][:,0:301]
                tab_ok["FLUX_H_TOTAL_UNIF"]=tab_ok["FLUX_H_APER"]

                #names spe
                
                # add columns spe 
                #if add_spe:
                    #names_spe = tab_s.colnames[1:-1]
                    #for name in names_spe:
                    #    tab_ok.add_column(tab_s[sgal][name],copy=True)
                    
                # Selection of MER  galaxies:
                threshold=0.8
                filt = filter_gals(tab_ok, mask_files=mask_files, threshold=threshold)
                tab_ok = tab_ok[filt]
            
                if mask_files:
                    #area_pix,area_weighted, relative_area = make_unif_mask(mask_files, subdir_results, threshold)
                    dummy = make_unif_mask(mask_files, subdir_results, threshold)
                    print("mask done")
                else:
                    area_pix=0
                    relative_area=0
                    area_weighted=0

            #with open(stat_file, 'a') as statfile:
            #    statfile.write(f"{tile_nb} \t {np.sum(filt)} \t {len(filt)} \t {area_pix:.4f} \t {area_weighted:.4f} \t {relative_area:.2f} \n")

            if add_spe:    
                nfile = 'EUC_MER_PHZ_SPE_CAT'+f"_tresh{threshold}_"+ fm.split('EUC_MER_FINAL-CAT')[1]
            else:
                nfile = 'EUC_MER_PHZ_CAT'+f"_tresh{threshold}_"+ fm.split('EUC_MER_FINAL-CAT')[1]
            #tab_ok.write('dataMERPHZ/'+nfile,format='fits')
            """
            tree = ET.parse(mer_xmls[midx])
            root = tree.getroot()
            newfile = root.find('Data/DataStorage/DataContainer/FileName')
            newfile.text = nfile
            #nxmlfile =  dir_files+subdir_results+'MERPHZCatalog'+mer_xmls[midx].split('DpdMerFinalCatalog')[1]
            #tree.write(nxmlfile)
            """
            
            colnames = np.intersect1d(tab_ok.colnames,cluster_fields)
            tab_cluster = tab_ok[colnames.tolist()]
            #pdb.set_trace()
            # nfile = 'EUC_CL_CAT'+ fm.split('EUC_MER_FINAL-CAT')[1]
            tab_cluster.write(dir_files+subdir_results+nfile,format='fits',overwrite=True)
            print(dir_files+subdir_results+nfile)
            tree = ET.parse(mer_xmls[midx])
            root = tree.getroot()
            newfile = root.find('Data/DataStorage/DataContainer/FileName')
            newfile.text = nfile
            if add_spe:
                nxmlfile = dir_files+subdir_results+"MERPHZSPECatalog"+mer_xmls[midx].split('DpdMerFinalCatalog')[1]
            else:
                nxmlfile = dir_files+subdir_results+"MERPHZCatalog"+mer_xmls[midx].split('DpdMerFinalCatalog')[1]
            tree.write(nxmlfile)
            print("Writting "+nxmlfile)
            

####################################################"

dir_files='/sps/euclid/OU-LE3/VMPZ-ID/newF006/sandbox/'
#dir_files='/sps/euclid/OU-LE3/VMPZ-ID/newF006/'
dir_files='/sps/euclid/OU-LE3/VMPZ-ID/REGREPROC1_R2/'
dir_files='/sps/euclid/OU-LE3/VMPZ-ID/Q1/'
phz_str = 'DpdPhzPfOutputCatalog*xml'
mer_str = 'DpdMerFinalCatalog-MER*xml'
spe_str = 'DpdSpePfOutputCatalog*xml'
subdir_results='CatRed_CL/'

#mask_fits='EUC_LE3_VMPZ-ID_HPCOVERAGE-%i_20240409T210514.345836Z_0.10.fits' # %i is the tile_index 
mask_str='DpdHealpixCoverageVMPZ*.xml'
import pdb

if __name__=="__main__":
    #print(dir_files)
    main_merge(dir_files, phz_str, mer_str, spe_str, mask_str, subdir_results=subdir_results,add_spe=False)

