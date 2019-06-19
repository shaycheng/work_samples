#!/usr/bin/env python 
'''
@author scheng3 
shuiting.cheng@ucsf.edu
'''
import os
import sys
import glob
from subprocess import call
from dipy.align.reslice import reslice
import nibabel as nib
from scipy import ndimage
from scipy.ndimage import label
import numpy as np
from dipy.io.image import load_nifti, save_nifti


mse_ID_list = sys.argv[1]


def dil_fast(CZ4, DSC):
    # check DSCv1, DSCv1 *(Fpv1 + DSC) = DSCF1, DSCv2 *(Fpv2 + DSC) = DSCF2
    DSCv1 = DSC[:-7] + '_v1.nii.gz'
    DSCv2 = DSC[:-7] + '_v2.nii.gz'
    DSCv3 = DSC[:-7] + '_v3.nii.gz'
    call(['fslmaths', DSC, '-kernel', '2D', '-dilM', DSCv1])
    call(['fslmaths', DSCv1, '-kernel', '2D', '-dilM', DSCv2])
    call(['fslmaths', DSCv2, '-kernel', '2D', '-dilM', DSCv3])
    CZ4v3 = CZ4[:-7] + '_v3.nii.gz'
    call(['fslmaths', CZ4, '-mas', DSCv3, CZ4v3])
    CZ4v3n3 = CZ4v3[:-7] + '_n3.nii.gz' 
    call(['fast', '-n', '3', '-t', '2', '-o', CZ4v3n3, CZ4v3])   
    # get DSCv1 *(Fpv1 + DSC) => DSCF1
    DSCF1 = DSC[:-7] + '_DSCF1.nii.gz'
    Fpv1 = glob.glob(os.path.join(path, '*v3_n3_pve_1.nii.gz'))
    call(['fslmaths', DSC, '-add', Fpv1[0], DSCF1])
    call(['fslmaths', DSCv1, '-mul', DSCF1, DSCF1])
    # get DSCv2 *(Fpv2 + DSC) => DSCF2
    DSCF2 = DSC[:-7] + '_DSCF2.nii.gz'
    Fpv2 = glob.glob(os.path.join(path, '*v3_n3_pve_2.nii.gz'))
    call(['fslmaths', DSC, '-add', Fpv2[0], DSCF2])
    call(['fslmaths', DSCv2, '-mul', DSCF2, DSCF2])
    call(['rm', CZ4v3])
    return DSCv1, DSCF1, DSCF2


def vessel_pop(DSC):
    # load DSC and get numbers of label
    cord_data, cord_affine = load_nifti(DSC)
    label_array,label_num=label(cord_data)
    print('cluster numbers from DSC are {}'.format(label_num)) 
    while label_num > 1:
         # while label number > 1, save labeled seg and get DSC_vessel, mask it to CZ4 and then sub such vessel volume from CZ4
        label_mask = DSC[:-7]+'_label.nii.gz'
        save_nifti(label_mask, label_array, cord_affine)    
        vessel = DSC[:-7] + '_vessel.nii.gz'
        call(['fslmaths', label_mask, '-thr', '2', vessel])
        CZ4_vessel = CZ4[:-7] + '_vessel.nii.gz'
        call(['fslmaths', CZ4, '-mas', vessel, CZ4_vessel])
        CZ4_no_vessel = CZ4[:-7] + '_no_vessel.nii.gz'
        call(['fslmaths', CZ4, '-sub', CZ4_vessel, CZ4_no_vessel])
        # run deepseg again based on CZ4_no_vessel
        qc = os.path.join(path, 'qc')
        call(['sct_deepseg_sc', '-i', CZ4_no_vessel, '-c', 't2', '-ofolder', path, '-qc', qc])
        call(['rm', '-r', qc])
        re_DSC = CZ4_no_vessel[:-7]+'_seg.nii.gz'
        data, affine = load_nifti(re_DSC)
        label_array,label_num=label(data)
        DSC = re_DSC
        print(DSC)
        print('final label_num is {}'.format(label_num))
        call(['rm', label_mask, vessel, CZ4_vessel])
    return DSC


def call_deepseg(CZ4):
    # extend CZ4 40 slices = 20 mm
    z_0 = CZ4[:-7]+'_z0.nii.gz'
    call(['fslroi', CZ4, z_0, '0', '-1', '0', '-1', '0', '1'])
    z_10 = CZ4[:-7] + '_z10.nii.gz'
    call(['fslmerge', '-z', z_10, z_0, z_0, z_0, z_0, z_0, z_0, z_0, z_0, z_0, z_0])
    z_40 = CZ4[:-7] + '_z40.nii.gz'
    call(['fslmerge', '-z', z_40, z_10, z_10, z_10, z_10])
    CZ4_extend = CZ4[:-7] + '_extend.nii.gz'
    call(['fslmerge', '-z', CZ4_extend, z_40, CZ4])
    call(['rm', z_0, z_10, z_40])
    # run deepseg on CZ4_extend
    qc = os.path.join(path, 'qc')
    call(['sct_deepseg_sc', '-i', CZ4_extend, '-c', 't2','-ofolder', path, '-qc', qc])
    call(['rm', '-r', qc])
    # abandon the appended 40 slices from deepseg and delete CZ4_extend & DSC_extend images
    DSC_extend = CZ4_extend[:-7]+'_seg.nii.gz'
    DSC = CZ4[:-7]+'_deepseg.nii.gz'
    call(['fslroi', DSC_extend, DSC, '0', '-1', '0', '-1', '40', '96'])
    call(['rm', DSC_extend, CZ4_extend])
    return DSC


def get_CZ4(path, T2W4):
    # call sct_get_centerline
    call(['sct_get_centerline', '-i', T2W4, '-c', 't2', '-ofolder', path])
    centerline = glob.glob(os.path.join(path, '*centerline_optic.nii.gz'))
    centerline = centerline[0]
    # locate cord center point (x, y, 0) of the bottom slice for crop 
    centerline_data, centerline_affine = load_nifti(centerline)
    center_point = np.nonzero(centerline_data[:, :, 0])
    x_center = np.array(center_point[0])
    x_left = x_center[0] - 20
    x_right = x_center[0] + 20 
    y_center = np.array(center_point[1])
    y_left = y_center[0] - 20
    y_right = y_center[0] + 20
    # crop T2W4 image from (x, y, 0)
    img = nib.load(T2W4)
    imgdata = img.get_data()
    imgaff = img.affine 
    img_crop = imgdata[x_left:x_right, y_left:y_right, 0:16]  
    # reslice T2W4_crop image to 0.5 mm3
    zooms = img.header.get_zooms()[:3]
    new_zooms = [0.5, 0.5, 0.5]
    newdata, newaff = reslice(img_crop, imgaff, zooms, new_zooms) 
    new_img_nii = nib.Nifti1Image(newdata, newaff)
    CZ = T2W4[:-7]+'_CZ.nii.gz'
    nib.save(new_img_nii, CZ)
    # call N4BiasFieldCorrection for CZ to get CZ4
    CZ4 = T2W4[:-7]+'_CZ4.nii.gz'
    call(['N4BiasFieldCorrection', '-i', CZ, '-o', CZ4])
    call(['rm', CZ])
    return CZ4


def create_T2W4(henryID):
    # copy T2W from PBR to local folder
    PBR_path = '/data/{}/PBR/subjects/{}/nii'.format(henryID, mseID)
    T2_name = ('t2_tse_tra_p2.nii.gz', 'T2_axial_tse_p3_3ml.nii.gz', 'T2_AX.nii.gz', 'T2_AX____60_SLICES.nii.gz', 'T2_AX____60_SLICES-000.nii.gz')     
    call(['mkdir', '/data/henry2/scheng/PROJECT/OPERA_t1t2pd/OPERA_t2/subjects/{}'.format(mseID)])
    path = '/data/henry2/scheng/PROJECT/OPERA_t1t2pd/OPERA_t2/subjects/{}'.format(mseID) 
    for each in os.listdir(PBR_path):
        if each.endswith(T2_name):
            t2 = os.path.join(PBR_path, each)
    call(['cp', t2, path])
    # call N4BiasFieldCorrection for T2W to get T2W4
    for each in os.listdir(path):
        if each.endswith(T2_name):
            T2W = os.path.join(path, each)
            T2W4 = T2W[:-7]+'_N4.nii.gz'
            call(['N4BiasFieldCorrection', '-i', T2W, '-o', T2W4])
    return path, T2W4


if __name__ == '__main__':
    i = 1
    with open(mse_ID_list, 'r') as M:
        for line in M:
            mseID = line[:-1]
            print(i)
            print(mseID)  
            num = int(''.join(filter(str.isdigit, mseID)))
            if num < 4700:   
                try:       
                    ### N4 bias correct T2W => T2W4   
                    path, T2W4 = create_T2W4('henry7')
                    ### Centerline_Crop-Zoom-N4 bias correction => CZ4
                    CZ4 = get_CZ4(path, T2W4)
                    ### Run DeepSeg-T2 on CZ4_extend to get cord => DSC
                    DSC = call_deepseg(CZ4)
                    ### make sure the cord_deepseg only contain cord seg, no vessel => DSC
                    DSC = vessel_pop(DSC)
                    ### dilation and run fast => DSCv1, DSCF1, DSCF2
                    DSCv1, DSCF1, DSCF2 = dil_fast(CZ4, DSC)
                except:
                    print('check {}'.format(mseID))
            else:             
                try:
                    ### N4 bias correct T2W => T2W4   
                    path, T2W4 = create_T2W4('henry11')
                    ### Centerline_Crop-Zoom-N4 bias correction => CZ4
                    CZ4 = get_CZ4(path, T2W4)
                    ### Run DeepSeg-T2 on CZ4_extend to get cord => DSC
                    DSC = call_deepseg(CZ4)
                    ### make sure the cord_deepseg only contain cord seg, no vessel => DSC
                    DSC = vessel_pop(DSC)
                    ### dilation and run fast => DSCv1, DSCF1, DSCF2
                    DSCv1, DSCF1, DSCF2 = dil_fast(CZ4, DSC)
                except:
                    print('check {}'.format(mseID))
            i = i + 1 
