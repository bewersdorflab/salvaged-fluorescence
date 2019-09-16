# SF_color_assignment
This package assigns single-molecule identities based on the salvaged fluorescence approach.

# Requirements
  - Microsoft Windows 7 or newer, 64-bit
  - Matlab R2017b or newer  
    - Curve Fitting Toolbox
    - Optimization Toolbox
  
# How to run
example1.m: this example shows how to estimate the photon number in the salvaged fluorescence channel.

example2.m: this example shows how to perform color assignment based on the photon numbers in the conventional and salvaged fluorescence channel.

# Optional step
Conventional_fitting.m: this script shows how to perform single-molecule fitting in the conventional fluorescence channel and estimate the positions and photon numbers. It saves the fitting result as "data_example1.mat", which will be loaded in example1.m. If the fitting result does not exist, it will run this script automatically to generate the data.

The optional step uses a GPU single-molecule fitter from "fit3Dcspline"(https://github.com/jries/fit3Dcspline; https://www.nature.com/articles/nmeth.4661)

The GPU fitter requires:
  - Microsoft Windows 7 or newer, 64-bit
  - CUDA capable graphics card, minimum Compute Capability 3.0
  - CUDA 8 compatible graphics driver (for GeForce products 378.66 or later)

# Data required
An example two-color dataset of a COS-7 cell imaged with the salvaged fluorescence approach. 

Download the dataset from: https://www.dropbox.com/sh/qw18jf3qcma1pbj/AACVAfhZ0WxM87vW-5VgrvIaa?dl=0
  - Raw images (3000 frames):  
    - Conventional_fluorescence.tif: images in the conventional fluorescence channel, acuiqred by a sCMOS camara
    - Salvaged_fluorescence.tif: images in the salvaged fluorescence channel, acuiqred by an EMCCD camara
  - Fitting results:
    - data_example1.mat: this is the fitting result of "conventional fluorescence.tif", generated by running "Conventional_fitting.m"
    - data_example2.mat: this is the analyzed result (provided by the authors) of the full dataset (150,000 frames) that used in example2.m to perform color assignment and reconstruction 

Additonal information of the dataset
  - Labeling: 
    - Channel 1: overexpressed GFP-Sec61b labeled with anti-GFP primary antibody and CF660C conjugated secondary antibody
    - Chennel 2: anti-a-tubulin primary antibody and AF647 conjugated secondary antibody
  - Imaging conditions: 
    - 200 fps 
    - 642 nm laser at 15 kW/cm2 
    - only the lower objective was used to collect the fluorescence

# Contact
For any questions / comments about this software, please contact [Bewersdorf Lab](http://www.bewersdorflab.org/).

The package is licenced under the [GNU GPL](https://www.gnu.org/licenses/). 

# Copyright and Software License
Copyright (c) 2019 Bewersdorf Lab, Yale Univeristy School of Medcine, USA.

# How to cite SF_color_assignment
If you use SF_color_assignment to process your data, please, cite our [paper](https://www.biorxiv.org/content/10.1101/613174v1):
  * Zhang, Y. et al. Nanoscale subcellular architecture revealed by multicolor 3D salvaged fluorescence imaging. bioRxiv, 613174 (2019).
 
