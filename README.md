# HPP

A multi-modal approach that combines structural and functional connectivity, termed hybrid pallidum parcellation (HPP), to parcellate the GP into anterior (aGP), middle (mGP), and posterior (pGP) subregions



## Project Overview

This MATLAB-based project provides a complete pipeline for processing resting-state fMRI (rs-fMRI) and diffusion tensor imaging (DTI) data to analyze subcortical brain structures, such as the Globus Pallidus (GP). The pipeline includes data preprocessing, functional connectivity (FC) computation, structural connectivity (SC) integration, density-based clustering (SNNDPC) for parcellation, generation of random parcellations as null models, and calculation of multimodal connectivity homogeneity (using PCA on combined SC+FC profiles).



The primary goal is to unveil topographic organization in subcortical ROIs by combining FC gradients with SC, perform clustering to identify parcels, and assess their homogeneity against random null models for statistical significance. This is inspired by methods like those in Tian et al. (2020) for subcortical mapping and is suitable for neuroimaging research on brain networks, e.g., in Parkinson's disease or cognitive studies.



Key features:

- Handles HCP-style CIFTI data for fMRI.

- Integrates DTI-based SC matrices.

- Generates null models for hypothesis testing.

- Outputs include clustered brain maps (NIfTI) and homogeneity scores.



## Dependencies

The pipeline requires the following software and toolboxes. Install them before running the scripts:



- **MATLAB R2023b** (or compatible version; tested on R2023b for optimal performance with built-in functions like `pca` and `pdist`).

- **Connectome Workbench**: For `wb_command` (used in spatial smoothing). Download from [Human Connectome Project](https://www.humanconnectome.org/software/connectome-workbench) and add to your system's PATH.

- **cifti-matlab-master**: For `cifti_read` (loading CIFTI files). Clone from [GitHub](https://github.com/Washington-University/cifti-matlab) and add to MATLAB path.

- **DPABI Toolbox**: For `y_IdealFilter` (band-pass filtering). Download from [RFMRI.org](http://rfmri.org/DPABI) and add to path.

- **Tian et al. (2020) Subcortex Functions**: For `read` (NIfTI loading) and `mat2nii` (matrix to NIfTI conversion). Clone the functions folder from [GitHub](https://github.com/yetianmed/subcortex/tree/master/functions) and add to MATLAB path. Reference: Tian et al. (2020) "Topographic organization of the human subcortex unveiled with functional connectivity gradients" (Nature Neuroscience).



Additional notes:

- MATLAB's Statistics and Machine Learning Toolbox is required for `pdist`, `pca`, `corr`, etc. (included in standard installations).

- Java support is needed for `java.util.LinkedList` in SNNDPC (built-in to MATLAB).

- No internet access is needed during runtime, but ensure all toolboxes are downloaded and paths added via `addpath`.



## Usage Guide

The project consists of three main scripts. Run them sequentially for a full analysis. Each script has a CONFIGURATION section for paths and parameters—customize as needed.



### 1. Preprocessing_SNNDPC.m

**Description**: Performs fMRI preprocessing (smoothing, filtering, normalization), extracts ROI time series, computes FC using PCA, aggregates group-level similarity matrices, combines with DTI data, and applies SNNDPC clustering to generate parcellations. Outputs include time series files, FC matrices, and clustered NIfTI maps (e.g., `GP_left3T.nii`).



### 2. Null_Model_Generation.m

**Description**: Generates 100 random spatial parcellations within left/right ROIs as a null model for statistical testing (e.g., compare homogeneity against random distributions).



### 3. Homogeneity.m

**Description**: Computes multimodal connectivity homogeneity for each parcel in the ROI across subjects, using PCA on combined SC+FC profiles. High scores indicate consistent connectivity within parcels.



**Reference Version**: `Homogeneity_me.m`, `Preprocessing_SNNDPC_me.m`, `Null_Model_Generation_me.m`includes example user-specific paths (e.g., `E:\research\GP-sub\...`). Use it as a template to replace placeholders in `Homogeneity.m`,  `Preprocessing_SNNDPC.m`, `Null_Model_Generation.m`, but the clean version is for general use.



## Notes and Troubleshooting

- **Path Replacement**: All scripts have placeholder paths—replace with your local directories (e.g., based on `Homogeneity_me.m` examples). Use absolute paths for reliability.

- **ROI Assumptions**: Scripts assume bilateral ROIs (left/right) with 3 parcels each; adjust for your atlas.

- **Performance**: For large cohorts, add `parfor` to loops; monitor memory for big matrices.

- **Errors**: If "undefined function" (e.g., read), ensure toolboxes are added to path. Test with sample data.



For issues, check MATLAB console logs or tool-specific docs.

