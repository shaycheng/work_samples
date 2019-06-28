## Pipelines for Processing Brain and Spinal Cord MRI data 
* dMRI pipeline * https://github.com/shaycheng/work_samples/blob/master/tutorial_dMRI_pipeline/dMRI_pipeline.ipynb
 ![brain_bundle](/tutorial_dMRI_pipeline/cst_track.png)
 
 
    The first one is to get matric values from specific region of interests in spinal cord C2-C3 sections from 
    multimodal images as features for later on statistical analysis. It involves registraion, transformation,
    crop and upsampling images and use deep learning technique to get total cord mask then use it to get image 
    ROI values. Check out script and result: *****C23_pipeline.py, C23_output.png*****
    
   
![C23_output](C23_output.png)




#### Required Software and Tools

Suppose you already have numpy, scipy, pandas, scikit-learn, seabron pip or conda installed

FSL                https://fsl.fmrib.ox.ac.uk/fsl/fslwiki

dipy               https://nipy.org/dipy/

spinalcord toolbox https://github.com/neuropoly/spinalcordtoolbox



