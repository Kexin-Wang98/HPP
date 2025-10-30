%% Create Spatially Constrained Random Parcellations (Null Model Generation)
%
% This script generates a set of random parcellations within predefined Regions 
% of Interest (ROIs), typically a left and a right hemisphere structure. The
% method used is a form of Voronoi tessellation based on Euclidean distance.
%
% The primary purpose of this script is to create a null model for statistical
% testing. By comparing a property of an empirically derived parcellation 
% (e.g., functional homogeneity) against the distribution of that property from
% these 100 random parcellations, one can assess its statistical significance.
%
% ## Pipeline
% 1.  **Load ROI Mask**: Reads a 3D volume (e.g., a NIfTI file) where different
%     integer values define the left and right ROIs.
% 2.  **Calculate Voxel Distances**: Computes the pairwise Euclidean distance
%     between all voxels within each ROI.
% 3.  **Iterative Parcellation**: For a specified number of iterations:
%     a. Randomly selects a set of "seed" voxels within each ROI.
%     b. Assigns every voxel in the ROI to the closest seed, creating parcels.
% 4.  **Store Results**: Saves the parcel assignments for each of the 100 iterations.
%
% ## Dependencies
% - MATLAB
% - A function to read 3D volume data (e.g., a NIfTI reader). The placeholder
%   `your_nifti_read_function` is used below.
%
%

%% ========================================================================
%  1. CONFIGURATION AND SETUP
% =========================================================================

% --- User-Defined Parameters ---

% Path to the NIfTI file containing the ROI mask.
% In this file, different integer labels should define the ROIs.
path_to_roi_mask = 'E:\research\GP-sub\7T_result\dtiAfmri\GP_random6.nii';

% Labels within the mask file that correspond to the left and right ROIs.
% For example, voxels with values 1, 2, or 3 belong to the left ROI.
labels_left_roi = [1, 2, 3];
labels_right_roi = [4, 5, 6];

% The desired number of parcels to create within EACH ROI.
num_parcels_per_roi = 3;

% The number of random parcellations to generate.
num_random_iterations = 100;

% The dimensions of the 3D volume (from the mask file).
% This is needed to convert linear voxel indices to 3D coordinates.
volume_dimensions = [113, 136, 113];

% --- End of Parameters ---


%% ========================================================================
%  2. LOAD DATA AND PREPARE COORDINATES
% =========================================================================

fprintf('Loading ROI mask and preparing voxel coordinates...\n');

% Load the ROI mask volume using a custom or toolbox function.
[~, roi_volume] = read(path_to_roi_mask);

% Get the linear indices of all voxels belonging to the left and right ROIs.
ind_left = find(ismember(roi_volume, labels_left_roi));
ind_right = find(ismember(roi_volume, labels_right_roi));

% Convert the linear indices of the ROI voxels into 3D [X, Y, Z] coordinates.
% This is essential for calculating spatial distances.
[coords_left(:,1), coords_left(:,2), coords_left(:,3)] = ind2sub(volume_dimensions, ind_left);
[coords_right(:,1), coords_right(:,2), coords_right(:,3)] = ind2sub(volume_dimensions, ind_right);


%% ========================================================================
%  3. PRE-CALCULATE PAIRWISE DISTANCE MATRICES
% =========================================================================
% To speed up the iterative process, we pre-calculate the Euclidean distance
% between every pair of voxels within each ROI.

fprintf('Calculating pairwise Euclidean distances for each ROI...\n');

% Calculate the distance matrix for the left ROI.
% 'pdist' returns a vector; 'squareform' converts it to a symmetric matrix.
dist_matrix_left = squareform(pdist(coords_left, 'euclidean'));

% Calculate the distance matrix for the right ROI.
dist_matrix_right = squareform(pdist(coords_right, 'euclidean'));


%% ========================================================================
%  4. GENERATE RANDOM PARCELLATIONS
% =========================================================================
% This is the main loop where the random parcellations are created.

fprintf('Generating %d random parcellations...\n', num_random_iterations);

% Pre-allocate matrices to store the results.
% Each column will represent one random parcellation.
parcel_assignments_left = zeros(size(coords_left, 1), num_random_iterations);
parcel_assignments_right = zeros(size(coords_right, 1), num_random_iterations);

for i = 1:num_random_iterations
    % --- Left ROI Parcellation ---
    % 1. Randomly select 'num_parcels_per_roi' voxels to act as parcel centers/seeds.
    center_indices_left = randperm(size(coords_left, 1), num_parcels_per_roi);
    
    % 2. Assign each voxel to the nearest center.
    % This is a fast way to perform a Voronoi-like tessellation. For each voxel
    % (row), 'min' finds the closest center (column) from the pre-computed
    % distance matrix and returns its index (i.e., the parcel label).
    [~, parcel_assignments_left(:, i)] = min(dist_matrix_left(:, center_indices_left), [], 2);
    
    % --- Right ROI Parcellation ---
    % 1. Repeat the process for the right ROI.
    center_indices_right = randperm(size(coords_right, 1), num_parcels_per_roi);
    
    % 2. Assign each voxel in the right ROI to its nearest random center.
    [~, parcel_assignments_right(:, i)] = min(dist_matrix_right(:, center_indices_right), [], 2);
    
    if mod(i, 10) == 0
        fprintf('  ... iteration %d of %d complete.\n', i, num_random_iterations);
    end
end

