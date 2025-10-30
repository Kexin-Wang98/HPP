%% Multimodal Connectivity Homogeneity Analysis (SC+FC)
%
% This script calculates the homogeneity of connectivity profiles for parcels
% within a specific Region of Interest (ROI), likely a subcortical structure 
% like the Globus Pallidus (GP). Homogeneity is defined as the variance 
% explained by the first principal component of a combined structural (SC) and 
% functional (FC) connectivity matrix.
%
% This analysis is performed for each parcel within the ROI and for each 
% subject in a dataset. High homogeneity suggests that the voxels within a 
% parcel share a similar pattern of brain-wide connectivity.
%
% ## Pipeline
% 1.  **Setup**: Defines paths, loads necessary brain masks (parcellations, ROIs),
%     and identifies voxel indices for each structure.
% 2.  **Subject Iteration**: Loops through each subject in the dataset.
% 3.  **Data Loading**: For each subject, loads their SC matrix (e.g., from DTI)
%     and FC time series data.
% 4.  **Parcel Iteration**: For each parcel within the ROI:
%     a. **Calculate SC Profile**: Computes the average structural connectivity
%        from the parcel's voxels to a set of cortical targets (e.g., Schaefer atlas).
%     b. **Calculate FC Profile**: Computes the functional connectivity (correlation)
%        between the parcel's time series and the time series of the same cortical targets.
%     c. **Combine and Analyze**: Concatenates the SC and FC profiles and runs a
%        Principal Component Analysis (PCA).
%     d. **Store Homogeneity**: The percentage of variance explained by the first
%        principal component is stored as the homogeneity score.
%
% ## Dependencies
% - MATLAB (Statistics and Machine Learning Toolbox for 'pca', 'corr')
% - Functions for reading neuroimaging data (e.g., NIfTI files). 
%
%

%% ========================================================================
%  1. CONFIGURATION AND INITIAL SETUP
% =========================================================================

% --- User-Defined Parameters ---

% Base directory for data and masks.
path_base_dir = 'E:\research\GP-sub\7T_result\';

% Paths to specific data types.
path_structural_data = 'E:\research\GP-sub\7T_result\new_result\result_nonline\nonline_diff_anal\';
path_functional_data = 'D:\function7Tdata_new\';

% Paths to atlas and mask files.
path_schaefer_atlas = 'E:\research\GP-sub\mask\Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1_6mm.nii';
path_yeo_networks = 'E:\research\GP-sub\mask\yeo17network_mask.nii';
path_roi_parcellation = 'E:\research\GP-sub\7T_result\dtiAfmri\GP_random6.nii'; % e.g., GP_random6.nii
path_roi_left_mask = 'E:\research\GP-sub\mask\GP_leftmask.nii';
path_roi_right_mask = 'E:\research\GP-sub\mask\GP_rightmask.nii';

% The number of parcels defined within EACH hemisphere of the ROI.
num_parcels_per_hemisphere = 3; 

% --- End of Parameters ---


%% ========================================================================
%  2. LOAD ATLASES AND DEFINE VOXEL INDICES
% =========================================================================
% This section loads all necessary masks and atlases and pre-calculates the
% linear indices for each region to speed up later processing.

fprintf('Loading atlases and defining voxel indices...\n');

% Load the cortical parcellation (e.g., Schaefer 400).
[~, schaefer_atlas] = read(path_schaefer_atlas);
% Load a network mask (e.g., Yeo 17) to constrain the atlas.
[~, yeo_networks] = read(path_yeo_networks);
% Apply the Yeo mask to the Schaefer atlas.
schaefer_atlas = schaefer_atlas .* yeo_networks; 
ind_cortical_voxels = find(yeo_networks ~= 0); % All voxels within the cortex mask

% Load the custom ROI parcellation (e.g., 6 parcels in the GP).
[~, roi_parcels] = read(path_roi_parcellation);
ind_roi_all_voxels = find(roi_parcels ~= 0); % All voxels within the entire ROI

% Load the binary masks for the left and right sides of the ROI.
[~, roi_left_mask] = read(path_roi_left_mask);
ind_roi_left_voxels = find(roi_left_mask ~= 0);

[~, roi_right_mask] = read(path_roi_right_mask);
ind_roi_right_voxels = find(roi_right_mask ~= 0);

% Get the list of subject directories.
File = dir(path_functional_data); 
FileNames = {File.name}';            
Length_Names = size(FileNames, 1);

% Pre-allocate the final results matrix.
total_parcels = 2 * num_parcels_per_hemisphere;
num_subjects = Length_Names - 2;
homogeneity_results = zeros(total_parcels, num_subjects);


%% ========================================================================
%  3. MAIN PROCESSING LOOP: COMPUTE HOMOGENEITY FOR EACH SUBJECT
% =========================================================================

fprintf('Starting homogeneity calculation for %d subjects...\n', num_subjects);

% Loop through all subjects (indices start at 3 to skip '.' and '..').
for i = 3:3
    subject_id = FileNames{i};
    fprintf('--- Processing Subject: %s (%d/%d) ---\n', subject_id, i-2, num_subjects);

    % --- Load Subject-Specific Data ---
    % Load structural connectivity matrix (ROI-to-cortex).
    sc_matrix_path = fullfile(path_structural_data, subject_id, 'fdt_matrix2over2.mat');
    load(sc_matrix_path, 'M2'); % Assumes the matrix is named M2

    % Load pre-processed functional time series.
    load(fullfile(path_functional_data, subject_id, 'GPL_TS.mat'), 'GPL_TS'); % Time serise of Left ROI
    load(fullfile(path_functional_data, subject_id, 'GPR_TS.mat'), 'GPR_TS'); % Time serise of Right ROI
    load(fullfile(path_functional_data, subject_id, 'net_ts_avg_400.mat'), 'net_ts_avg_400'); % Time serise of Cortical parcels in schaefer 400

    % --- Process Left Hemisphere Parcels ---
    for parcel_num = 1:num_parcels_per_hemisphere
        % Find the indices for voxels in the current parcel.
        ind_current_parcel_global = find(roi_parcels == parcel_num);
        % Map these global indices to their position within the full ROI mask.
        [~, ind_in_full_roi, ~] = intersect(ind_roi_all_voxels, ind_current_parcel_global);

        % -- Structural Connectivity Profile --
        sc_profile = [];
        % Loop through the first 200 cortical parcels (left hemisphere).
        for schaefer_parcel_num = 1:200
            % Find indices for the current cortical parcel.
            ind_schaefer_parcel_global = find(schaefer_atlas == schaefer_parcel_num);
            % Map to their position within the cortex mask.
            [~, ind_in_cortex_mask, ~] = intersect(ind_cortical_voxels, ind_schaefer_parcel_global);
            
            % Calculate the mean connectivity from the ROI parcel to the cortical parcel.
            mean_connectivity = mean(M2(ind_in_full_roi, ind_in_cortex_mask), 2);
            sc_profile(schaefer_parcel_num, :) = mean_connectivity';
        end

        % -- Functional Connectivity Profile --
        % Map global parcel indices to their position within the left ROI time series matrix.
        [~, ind_in_left_roi_ts, ~] = intersect(ind_roi_left_voxels, ind_current_parcel_global);
        % Correlate the cortical time series with the ROI parcel's time series.
        fc_profile = corr(net_ts_avg_400, GPL_TS(:, ind_in_left_roi_ts));
        
        % -- Homogeneity Calculation --
        combined_connectivity = [fc_profile; sc_profile];
        [~, ~, ~, ~, explained] = pca(combined_connectivity, 'Rows', 'complete');
        homogeneity_results(parcel_num, (i-2)) = explained(1); % Store variance explained by 1st PC
    end
    
    % --- Process Right Hemisphere Parcels ---
    for parcel_num = (num_parcels_per_hemisphere + 1):(2 * num_parcels_per_hemisphere)
        % Repeat the entire process for the right hemisphere parcels.
        ind_current_parcel_global = find(roi_parcels == parcel_num);
        [~, ind_in_full_roi, ~] = intersect(ind_roi_all_voxels, ind_current_parcel_global);

        % -- Structural Connectivity Profile --
        sc_profile = [];
        % Loop through the second 200 cortical parcels (right hemisphere).
        for schaefer_parcel_num = 201:400
            ind_schaefer_parcel_global = find(schaefer_atlas == schaefer_parcel_num);
            [~, ind_in_cortex_mask, ~] = intersect(ind_cortical_voxels, ind_schaefer_parcel_global);
            
            mean_connectivity = mean(M2(ind_in_full_roi, ind_in_cortex_mask), 2);
            sc_profile(schaefer_parcel_num-200, :) = mean_connectivity';
        end

        % -- Functional Connectivity Profile --
        [~, ind_in_right_roi_ts, ~] = intersect(ind_roi_right_voxels, ind_current_parcel_global);
        fc_profile = corr(net_ts_avg_400, GPR_TS(:, ind_in_right_roi_ts));
        
        % -- Homogeneity Calculation --
        combined_connectivity = [fc_profile; sc_profile];
        [~, ~, ~, ~, explained] = pca(combined_connectivity, 'Rows', 'complete');
        homogeneity_results(parcel_num, (i-2)) = explained(1);
    end
    
    % Clear subject-specific variables to free up memory.
    clear M2 GPR_TS GPL_TS net_ts_avg_400 fc_profile sc_profile;
end

fprintf('--- Analysis complete. Results are stored in the "homogeneity_results" matrix. ---\n');