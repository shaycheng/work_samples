#!/usr/bin/env python
'''
@author scheng3
shuiting.cheng@ucsf.edu
'''
import os 
import sys
import glob
from subprocess import Popen, PIPE
from subprocess import call
import csv 
import numpy as np
import nibabel as nib


mse_ID_list = sys.argv[1]


def get_median_values(image, mask):
    cmd = ["fslstats", image, "-k", mask, "-P", "50"]
    proc = Popen(cmd, stdout=PIPE)
    P_50 = [l.decode("utf-8").split() for l in proc.stdout.readlines()]
    P_50 = str(P_50)
    P_50 = P_50.replace("'","").replace("[[","").replace("]]","")
    #print(P_50)
    return P_50

def get_mean_values(image, mask):
    cmd = ["fslstats", image, "-k", mask, "-M"]
    proc = Popen(cmd, stdout=PIPE)
    MEAN = [l.decode("utf-8").split() for l in proc.stdout.readlines()]
    MEAN = str(MEAN)
    MEAN = MEAN.replace("'","").replace("[[","").replace("]]","")
    #print(MEAN)
    return MEAN

def write_to_csv():
    CSV = '/data/henry2/scheng/PROJECT/PMD/C23_ECTRIMS/C23_DSC_values.csv'
    img = nib.load(DSC)
    imgdata = img.get_data()
    vox_num = 0 
    for i in range(200):
        for j in range(200):
            if imgdata[i][j][0] == 1:
                vox_num = vox_num + 1
    zooms = img.header.get_zooms()[:3]
    vox_size = zooms[0] * zooms[1]
    TCA = vox_size * vox_num

    PSIR_DSC_median = get_median_values(PSIR, DSC)
    MT_DSC_median = get_median_values(MT, DSC)
    MTR_DSC_median = get_median_values(MTR, DSC)       
    FA_DSC_median = get_median_values(FA, DSC)
    MD_DSC_median = get_median_values(MD, DSC)
    AD_DSC_median = get_median_values(AD, DSC)
    RD_DSC_median = get_median_values(RD, DSC)

    PSIR_DSC_mean = get_mean_values(PSIR, DSC)
    MT_DSC_mean = get_mean_values(MT, DSC)
    MTR_DSC_mean = get_mean_values(MTR, DSC)   
    FA_DSC_mean = get_mean_values(FA, DSC)
    MD_DSC_mean = get_mean_values(MD, DSC)
    AD_DSC_mean = get_mean_values(AD, DSC)
    RD_DSC_mean = get_mean_values(RD, DSC)

    new_row = [mseID, TCA, PSIR_DSC_median, MT_DSC_median, MTR_DSC_median, FA_DSC_median, MD_DSC_median, AD_DSC_median, RD_DSC_median, mseID, PSIR_DSC_mean, MT_DSC_mean, MTR_DSC_mean, FA_DSC_mean, MD_DSC_mean, AD_DSC_mean, RD_DSC_mean]
    with open(CSV, 'a') as f:
        writer = csv.writer(f)
        writer.writerow(new_row)
        print(new_row)


if __name__ == '__main__':
    with open(mse_ID_list, 'r') as M:
        for line in M:
            mseID = line[:-1]
            print(mseID)
            try: 
                PMD = (glob.glob('/data/henry7/PBR/subjects/{}/PMD_C23'.format(mseID))+glob.glob('/data/henry11/PBR/subjects/{}/PMD_C23'.format(mseID)))[0]   
                PSIR = os.path.join(PMD, 'PSIR_CZ.nii.gz')
                MTR = os.path.join(PMD, 'MTR_reg_CZ.nii.gz')
                MT = os.path.join(PMD, 'MT_reg_CZ.nii.gz')
                DSC = os.path.join(PMD, 'PSIR_CZ_seg.nii.gz')
                FA = os.path.join(PMD, 'dti_FA_reg_CZ.nii.gz')
                MD = os.path.join(PMD, 'dti_MD_reg_CZ.nii.gz')
                RD = os.path.join(PMD, 'dti_RD_reg_CZ.nii.gz')
                AD = os.path.join(PMD, 'dti_AD_reg_CZ.nii.gz')
                write_to_csv()
            except:
                print('check {}'.format(mseID))

