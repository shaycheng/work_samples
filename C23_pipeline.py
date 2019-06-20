#!/usr/bin/env python 
'''
@author scheng3
shuiting.cheng@ucsf.edu
'''
import os
import sys
import glob
from subprocess import call
from subprocess import check_call
from dipy.align.reslice import reslice
import nibabel as nib
from scipy import ndimage
import numpy as np
from dipy.io.image import load_nifti, save_nifti


mse_ID_list = sys.argv[1]


def call_deepseg(img):
    # run deepseg on img
    qc = os.path.join(PMD, 'qc')
    call(['sct_deepseg_sc', '-i', img, '-c', 't1','-ofolder', PMD, '-qc', qc])
    call(['rm', '-r', qc])
    DSC = img[:-7]+'_seg.nii.gz'
    return DSC


def crop_zoom_file(image):
    centerline = glob.glob(os.path.join(PMD, 'PSIR_centerline_optic.nii.gz'))
    centerline = centerline[0]
    print(centerline)
    # locate cord center point (x, y, 0) of the bottom slice for crop 
    centerline_data, centerline_affine = load_nifti(centerline)
    center_point = np.nonzero(centerline_data[:, :, 0])
    x_center = np.array(center_point[0])
    x_left = x_center[0] - 25
    x_right = x_center[0] + 25
    y_center = np.array(center_point[1])
    y_left = y_center[0] - 25
    y_right = y_center[0] + 25
    print(x_left, y_left)
    # reslice crop image to 0.78/4 mm3
    img = nib.load(image)
    imgdata = img.get_data()
    imgaff = img.affine 
    img_crop = imgdata[x_left:x_right, y_left:y_right, :] 
    zooms = img.header.get_zooms()[:3]
    print(zooms)
    new_zooms = [zooms[0]/4, zooms[1]/4, zooms[2]]
    print(new_zooms)
    newdata, newaff = reslice(img_crop, imgaff, zooms, new_zooms) 
    new_img_nii = nib.Nifti1Image(newdata, newaff)
    CZ = image[:-7]+'_CZ.nii.gz'
    nib.save(new_img_nii, CZ)
    return CZ


def sct_register_multimodal():
    ## 1. register MTR to PSIR space and apply tranformation to MT
    # get correct slices for MT, MTR then use sct_deepseg_sc, sct_get_centerline to get centerline and cord_seg for multimodal registration
    call(['fslroi', MT, MT, '0', '-1', '0', '-1', '31', '2'])
    call(['fslroi', MTR, MTR, '0', '-1', '0', '-1', '31', '2'])
    PSIR = os.path.join(PMD, 'PSIR.nii.gz')
    call(['fslmerge', '-z', PSIR, PSIR_raw, PSIR_raw])
    call(['fslmerge', '-z', MTR, MTR, MTR])
    call(['fslmerge', '-z', MT, MT, MT])
    call(['rm', PSIR_raw])
    call(['sct_get_centerline', '-i', PSIR, '-c', 't1', '-ofolder', PMD])
    call(['sct_get_centerline', '-i', MTR, '-c', 't1', '-ofolder', PMD])
    qc = os.path.join(PMD, 'qc')
    call(['sct_deepseg_sc', '-i', PSIR, '-c', 't1', '-ofolder', PMD, '-qc', qc])
    call(['sct_deepseg_sc', '-i', MTR, '-c', 't1', '-ofolder', PMD, '-qc', qc])
    call(['rm', '-r', qc])           
    #register MTR to PSIR space
    MTR = os.path.join(PMD, 'MTR.nii.gz')
    MTR_centerline = os.path.join(PMD, 'MTR_centerline_optic.nii.gz')
    MTR_seg = os.path.join(PMD, 'MTR_seg.nii.gz')
    PSIR = os.path.join(PMD, 'PSIR.nii.gz')
    PSIR_centerline = os.path.join(PMD, 'PSIR_centerline_optic.nii.gz')
    PSIR_seg = os.path.join(PMD, 'PSIR_seg.nii.gz')
    call(['sct_register_multimodal', '-i', PSIR, '-iseg', PSIR_seg, '-ilabel', PSIR_centerline, '-d', MTR, '-dseg', MTR_seg, '-dlabel', MTR_centerline, '-ofolder', PMD, '-param', 'step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,slicewise=1,iter=3'])
    MTR_reg = os.path.join(PMD, 'MTR_reg.nii.gz')
    #apply warp_MTR2PSIR to MT 
    MT = os.path.join(PMD, 'MT.nii.gz')
    mt_warp = os.path.join(PMD, 'warp_MTR2PSIR.nii.gz')
    MT_reg = os.path.join(PMD, 'MT_reg.nii.gz')
    call(['sct_apply_transfo', '-i', MT, '-d', PSIR, '-w', mt_warp, '-o', MT_reg, '-x', 'nn']) 
    call(['fslmaths', MT_reg, '-thr', '0', '-uthr', '0.2', MT_reg])  
    warp2 = (glob.glob('/data/henry7/PBR/subjects/{}/PMD_C23warp_PSIR2MTR.nii.gz'.format(mseID))+glob.glob('/data/henry11/PBR/subjects/{}/PMD_C23warp_PSIR2MTR.nii.gz'.format(mseID)))[0] 
    call(['rm', warp2])

    ## 2. register dMRI data to PSIR and apply tranformation to dti metrics fa, md, rd, ad
    # compute mean dMRI 
    dmri_mean = os.path.join(PMD, 'dmri_mean.nii.gz')
    call(['sct_maths', '-i', dMRI_raw, '-mean', 't', '-o', dmri_mean])
    # get centerline from mean dMRI data
    call(['sct_get_centerline', '-i', dmri_mean, '-c', 'dwi', '-ofolder', PMD])
    # segment SC on mean dMRI data
    qc = os.path.join(PMD, 'qc')
    call(['sct_deepseg_sc', '-i', dmri_mean, '-c', 'dwi', '-ofolder', PMD, '-qc', qc])
    call(['rm', '-r', qc])
    # multimodal registration
    dmri_seg = os.path.join(PMD, 'dmri_mean_seg.nii.gz')
    dmri_centerline = os.path.join(PMD, 'dmri_mean_centerline_optic.nii.gz')
    call(['fslroi', dmri_mean, dmri_mean, '0', '-1', '0', '-1', '5', '1'])
    call(['fslmerge', '-z', dmri_mean, dmri_mean, dmri_mean])
    call(['fslroi', dmri_seg, dmri_seg, '0', '-1', '0', '-1', '5', '1'])
    call(['fslmerge', '-z', dmri_seg, dmri_seg, dmri_seg])
    call(['fslroi', dmri_centerline, dmri_centerline, '0', '-1', '0', '-1', '5', '1'])
    call(['fslmerge', '-z', dmri_centerline, dmri_centerline, dmri_centerline])
    call(['sct_register_multimodal', '-i', PSIR, '-iseg', PSIR_seg, '-ilabel', PSIR_centerline, '-d', dmri_mean, '-dseg', dmri_seg, '-dlabel', dmri_centerline, '-ofolder', PMD, '-param', 'step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,slicewise=1,iter=3'])
    # compute DTI metrics using dipy
    dti_warp = os.path.join(PMD, 'warp_dmri_mean2PSIR.nii.gz')
    call(['sct_dmri_compute_dti', '-i', dMRI_raw, '-bval', bval_raw, '-bvec', bvec_raw, '-o', os.path.join(PMD, 'dti_')])
    fa = os.path.join(PMD, 'dti_FA.nii.gz')
    call(['fslroi', fa, fa, '0', '-1', '0', '-1', '5', '1'])
    call(['fslmerge', '-z', fa, fa, fa])  
    fa_reg = fa[:-7]+'_reg.nii.gz'
    call(['sct_apply_transfo', '-i', fa, '-d', PSIR, '-w', dti_warp, '-o', fa_reg, '-x', 'nn'])

    md = os.path.join(PMD, 'dti_MD.nii.gz')
    call(['fslroi', md, md, '0', '-1', '0', '-1', '5', '1'])
    call(['fslmerge', '-z', md, md, md])  
    md_reg = md[:-7]+'_reg.nii.gz'
    call(['sct_apply_transfo', '-i', md, '-d', PSIR, '-w', dti_warp, '-o', md_reg, '-x', 'nn'])

    rd = os.path.join(PMD, 'dti_RD.nii.gz')
    call(['fslroi', rd, rd, '0', '-1', '0', '-1', '5', '1'])
    call(['fslmerge', '-z', rd, rd, rd])  
    rd_reg = rd[:-7]+'_reg.nii.gz'
    call(['sct_apply_transfo', '-i', rd, '-d', PSIR, '-w', dti_warp, '-o', rd_reg, '-x', 'nn'])

    ad = os.path.join(PMD, 'dti_AD.nii.gz')
    call(['fslroi', ad, ad, '0', '-1', '0', '-1', '5', '1'])
    call(['fslmerge', '-z', ad, ad, ad])  
    ad_reg = ad[:-7]+'_reg.nii.gz'
    call(['sct_apply_transfo', '-i', ad, '-d', PSIR, '-w', dti_warp, '-o', ad_reg, '-x', 'nn'])
    call(['rm', os.path.join(PBR, 'PMD_C23warp_PSIR2dmri_mean.nii.gz'), dMRI_raw, bvec_raw, bval_raw, fa, md, rd, ad])
    return PSIR, MTR_reg, MT_reg, fa_reg, md_reg, rd_reg, ad_reg



def create_mtr_spine(image_MT, image_1, image_2):
    print('ctreating MT, MTR now................')
    alpha_MT = (5*2*3.14159)/360
    alpha_1 = (5*2*3.14159)/360
    alpha_2 = (15*2*3.14159)/360
    TR_MT = 29.000
    TR_1 = 29.000
    TR_2 = 29.000

    image_MT_reg = image_MT[:-7]+'_reg.nii.gz'
    image_1_reg = image_1[:-7]+'_reg.nii.gz'
    call(['flirt', '-dof', '6', '-in', image_MT, '-ref', image_2, '-out', image_MT_reg])
    call(['flirt', '-dof', '6', '-in', image_1, '-ref', image_2, '-out', image_1_reg])

    img_MT = nib.load(image_MT_reg)
    img_1 = nib.load(image_1_reg)
    img_2 = nib.load(image_2)
    signal_MT = img_MT.get_data()
    signal_1 = img_1.get_data()
    signal_2 = img_2.get_data()

    np.seterr(divide='ignore', invalid="ignore")

    R_num1 = np.zeros(np.shape(signal_1))
    R_num1 = (alpha_1/(2*TR_1))*signal_1
    R_num2 = (alpha_2/(2*TR_2))*signal_2 
    R_num = R_num1-R_num2
    R_denom = (signal_2/alpha_2)-(signal_1/alpha_1)
    R = R_num/R_denom
    R[np.isinf(R)] = 0
    R = np.nan_to_num(R)

    delta_term1_1 = signal_1/signal_MT                                                   
    delta_term1_1[np.isinf(delta_term1_1)] = 0
    delta_term1_1 = np.nan_to_num(delta_term1_1)
    delta_term1 = delta_term1_1 - np.ones(np.shape(signal_1))
    delta_term2_1 = np.zeros(np.shape(signal_1))
    delta_term2_1[:] = (alpha_1*alpha_1)/2
    delta_term2_2 = TR_1*R
    delta_term2 = delta_term2_1 + TR_1*R
    delta = delta_term1*delta_term2  
    MT_img = nib.Nifti1Image(delta.astype('float32'), img_1.get_affine())
    nib.save(MT_img, os.path.join(PMD, 'MT.nii.gz'))
    MT = glob.glob(os.path.join(PMD, 'MT.nii.gz'))
    MT = MT[0]
    
    ratio = (signal_1 - signal_MT)/signal_1
    ratio[np.isinf(ratio)] = 0
    ratio[np.isnan(ratio)] = 0
    MTR_img = nib.Nifti1Image(ratio.astype('float32'), img_1.get_affine())
    MTR = nib.save(MTR_img, os.path.join(PMD, 'MTR.nii.gz'))
    MTR = glob.glob(os.path.join(PMD, 'MTR.nii.gz'))
    MTR = MTR[0]
    call(['rm', image_MT_reg, image_1_reg, MT120_raw, NONMT_raw, NONMTFA_raw])
    return MT, MTR


def copy_file(name):
    for each in os.listdir(nii_path):
        if each.endswith(name):
            files = os.path.join(nii_path, each)
            call(['cp', files, PMD])
            img = glob.glob(os.path.join(PMD, '*{}'.format(name)))
            img = img[0]
    return img


if __name__ == '__main__':
    with open(mse_ID_list, 'r') as M:
        for line in M:
            mseID = line[:-1]
            print(mseID)
            try: 
                PBR = (glob.glob('/data/henry7/PBR/subjects/{}'.format(mseID))+glob.glob('/data/henry11/PBR/subjects/{}'.format(mseID)))[0]
                nii_path = os.path.join(PBR, 'nii')
                call(['mkdir', os.path.join(PBR, 'PMD_C23')])
                PMD = os.path.join(PBR, 'PMD_C23')
                ### get PSIR, MTR, DIFFUSION raw images ready to PMD_C23 folder
                PSIR_raw = copy_file('-C2_3_2Fl_seg_psir_TI_PSIR.nii.gz')
                MT120_raw = copy_file('fl3d_TRA_p2_iso_120mm.nii.gz')
                NONMT_raw = copy_file('fl3d_TRA_p2_iso_120mm_NONMT.nii.gz')
                NONMTFA_raw = copy_file('fl3d_TRA_p2_iso_120mm_FA15_NONMT.nii.gz')
                dMRI_raw = copy_file('diff_tra_4mm_ZOOMit_dti_6dir_DFC-000.nii.gz')           
                bvec_raw = copy_file('diff_tra_4mm_ZOOMit_dti_6dir_DFC-002.bvec')            
                bval_raw = copy_file('diff_tra_4mm_ZOOMit_dti_6dir_DFC-001.bval')
                ### create MTR and MT 
                MT, MTR = create_mtr_spine(MT120_raw, NONMT_raw, NONMTFA_raw)
                print(MT, '\n', MTR)
                ### multimodal images registration to PSIR space
                PSIR, MTR_reg, MT_reg, fa_reg, md_reg, rd_reg, ad_reg = sct_register_multimodal()
                ### crop and zoom images 
                PSIR_CZ = crop_zoom_file(PSIR)
                crop_zoom_file(MTR_reg)
                crop_zoom_file(MT_reg)               
                crop_zoom_file(fa_reg)
                crop_zoom_file(md_reg)
                crop_zoom_file(rd_reg)
                crop_zoom_file(ad_reg)
                ### run deepseg on PSIR_CZ to get total cord mask
                DSC = call_deepseg(PSIR_CZ)
                call(['rm', fa_reg, md_reg, rd_reg, ad_reg])  
            except:
                print('check {}'.format(mseID))
    
