# **LDCT Simulation**  

**Low-Dose CT (LDCT) Simulation and Analysis**  

This repository contains the code used for a project focused on simulating and analyzing Low-Dose Computed Tomography (LDCT) data.

## **Dataset Sources**  
- **Primary Datasets:** 
[Colorectal-Liver-Metastases](https://www.cancerimagingarchive.net/collection/colorectal-liver-metastases/)
[LDCT-and-Projection-Data (The Cancer Imaging Archive)](https://www.cancerimagingarchive.net/collection/ldct-and-projection-data/)  
- **Supplementary Data:** Custom-generated phantom data  

## **Get Started**

This repository contains a folder called **`algorithm`**, which includes two Python files:

1. **`sim_functions.py`**:
   This file contains the main algorithm used to simulate image noise. It includes the function *`apply_ld_sinogram_full`*, as well as the function *`metrics_simulation_phan`*, which was used to generate the low-dose versions of the phantoms and to compute the SNR and CNR metrics.

2. **`utils_functions.py`**:
   This file contains utility functions that facilitate result visualization, metric computation, and other supportive tasks related to the main algorithm.

Outside the folder, the file **`save_colorectal_ldct.py`** is provided. This script was used to save the low-dose versions of the CRLM dataset at half dose, quarter dose, and one-eighth of the original dose.

There is also the file **`try_the_algorithm.ipynb`**, a Jupyter notebook that contains a step-by-step, explicit description of the algorithm, where its functionality can be tested using a CT slice and the results can be visualized.

Finally, there is the file **`solution_noise_sim.ipynb`**, which was used to save the simulated versions of the phantoms obtained at Dose 3.

To use the code provided in this repository, it is recommended to use a virtual environment and install all the required libraries, which are listed in the **`requirements.txt`** file.

## **wip [Work in progress]**

The file **`calculate_seg_tumor.ipynb`** corresponds to an ongoing work whose objective is to use the TotalSegmentator model to compare tumor segmentations across images acquired at different dose levels, generated using the simulation implemented in this repository.
However, during the implementation, issues related to the correct usage of the model were encountered, and therefore it has not yet been possible to obtain reliable results or complete this analysis.


## **References**   
Zeng, D., Huang, J., Bian, Z., Niu, S., Zhang, H., Feng, Q., Liang, Z., & Ma, J. (2015).  
*A Simple Low-Dose X-ray CT Simulation from High-Dose Scan.*  
**IEEE Transactions on Nuclear Science**, *62*(5), 2226–2233.    
