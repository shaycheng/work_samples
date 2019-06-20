## Pipelines for Processing Brain and Spinal Cord MRI data 
There are two pipeline samples here 
    
    The first one is to get matric values from specific region of interests in spinal cord C2-C3 sections from multimodal 
    images as  features for later on statistical analysis. It involves registraion, transformation, crop and upsampling 
    images and use deep learning technique to get total cord mask to get ROI values. Check out pipeline script and 
    output:  C23_pipeline, C23_output


for spinal cord total area measurement


#### Required Software and Tools


Suppose you already have numpy, scipy, pandas, scikit-learn, seabron pip or conda installed

FSL                https://fsl.fmrib.ox.ac.uk/fsl/fslwiki

dipy               https://nipy.org/dipy/

spinalcord toolbox https://github.com/neuropoly/spinalcordtoolbox



