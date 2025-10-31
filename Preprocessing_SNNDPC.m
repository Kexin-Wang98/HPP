%% fMRI Preprocessing and Functional Connectivity Pipeline
%
% This script performs a complete fMRI analysis pipeline for a group of subjects.
% The main steps are:
%   1. Spatial smoothing of CIFTI functional data.
%   2. Temporal preprocessing: discarding initial frames, band-pass filtering, 
%      and normalization.
%   3. Functional Connectivity (FC) analysis: calculating the connectivity between
%      pre-defined Regions of Interest (ROIs) and whole-brain principal components.
%   4. Aggregation of results across all subjects.
%
% Dependencies:
%   - MATLAB
%   - Connectome Workbench (must be in your system's PATH for 'wb_command').
%   - cifti-matlab-master (e.g., cifti_read).
%   - Custom filtering functions from DPABI (e.g., y_IdealFilter).
%   - Functions from Tian et al. (2020) 'Topographic organization of the human subcortex unveiled with functional connectivity gradients' 
%     (https://github.com/yetianmed/subcortex/tree/master/functions)
%
%
% Before running, please ensure all paths in the "CONFIGURATION" section are
% correctly set to your local directories.
%

%% ========================================================================
%  PART 1: SPATIAL SMOOTHING OF CIFTI DATA
% =========================================================================
% This section iterates through each subject's fMRI runs and applies spatial
% smoothing using Connectome Workbench's wb_command.
% -------------------------------------------------------------------------

% --- CONFIGURATION ---
% Define the root directories for your data and atlases.
% It is highly recommended to use relative paths or define all paths here.
path_base_input = '/path/to/your/project/directory';
path_subject_list = fullfile(path_base_input, 'subject_list_dir'); % Dir containing subject folders
path_raw_fmri_data = fullfile(path_base_input, 'raw_fmri_data'); % Parent dir for raw CIFTI files
path_smoothed_fmri = fullfile(path_base_input, 'derivatives', 'smoothed_fmri'); % Output for smoothed data
path_surface_atlases = fullfile(path_base_input, 'atlases', 'surface'); % Dir with .surf.gii files

% --- SETUP ---
File = dir(path_subject_list);
FileNames = {File.name}';
Length_Names = size(FileNames, 1);

% Define the smoothing kernel (Full Width at Half Maximum) in mm.
fwhm_kernel = 4; % for 7T data

% Start a parallel pool to speed up processing. Adjust the number of workers as needed.
parpool("local", 4);

% Loop through all subject folders (starting from 3 to skip '.' and '..').
parfor i = 3:Length_Names
    subject_id = FileNames{i};
    fprintf('--- Starting smoothing for subject: %s ---\n', subject_id);
    
    % Create an output directory for the current subject if it doesn't exist.
    out_dir = fullfile(path_smoothed_fmri, subject_id);
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % Define the four fMRI runs to be processed.
    run_names = { ...
        'rfMRI_REST1_7T_PA_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii', ...
        'rfMRI_REST2_7T_AP_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii', ...
        'rfMRI_REST3_7T_PA_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii', ...
        'rfMRI_REST4_7T_AP_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii'  ...
    };
    
    % Process each of the four runs for the subject.
    for j = 1:numel(run_names)
        input_file = fullfile(path_raw_fmri_data, subject_id, run_names{j});
        output_file = fullfile(out_dir, strrep(run_names{j}, '.dtseries.nii', '_smooth.dtseries.nii'));
        
        % Check if the raw input file exists before attempting to smooth.
        if exist(input_file, 'file')
            % Construct the wb_command for CIFTI smoothing.
            % This command smooths the data along the cortical surface, which is the
            % standard and recommended procedure for HCP-style data.
            left_surf = fullfile(path_surface_atlases, 'L.sphere.59k_fs_LR.surf.gii');
            right_surf = fullfile(path_surface_atlases, 'R.sphere.59k_fs_LR.surf.gii');
            
            cmd = sprintf(['wb_command -cifti-smoothing "%s" %d %d COLUMN "%s" ' ...
                '-left-surface "%s" -right-surface "%s"'], ...
                input_file, fwhm_kernel, fwhm_kernel, output_file, left_surf, right_surf);
            
            % Execute the command.
            system(cmd);
        else
            fprintf('WARNING: Input file not found for subject %s. Skipping: %s\n', subject_id, run_names{j});
        end
    end
    fprintf('Finished smoothing for subject: %s (%d/%d)\n', subject_id, i-2, Length_Names-2);
end


%% ========================================================================
%  PART 2: TEMPORAL PREPROCESSING
% =========================================================================
% This section loads the smoothed data and performs several temporal
% preprocessing steps:
%   - Discarding the first 15 volumes to avoid T1 saturation effects.
%   - Applying a band-pass filter to isolate the BOLD signal frequency.
%   - Normalizing the time series (z-scoring).
%   - Extracting time series for specific Regions of Interest (ROIs).
%   - Concatenating data from all four runs.
% -------------------------------------------------------------------------

% --- CONFIGURATION ---
path_processed_fmri = fullfile(path_base_input, 'derivatives', 'processed_fmri_timeseries'); % Final output dir
if ~exist(path_processed_fmri, 'dir'), mkdir(path_processed_fmri); end

% Define the smoothed file names (suffix added in Part 1).
run_templates = { ...
    '\rfMRI_REST1_7T_PA_Atlas_1.6mm_MSMAll_hp2000_clean_smooth.dtseries.nii', ...
    '\rfMRI_REST2_7T_AP_Atlas_1.6mm_MSMAll_hp2000_clean_smooth.dtseries.nii', ...
    '\rfMRI_REST3_7T_PA_Atlas_1.6mm_MSMAll_hp2000_clean_smooth.dtseries.nii', ...
    '\rfMRI_REST4_7T_AP_Atlas_1.6mm_MSMAll_hp2000_clean_smooth.dtseries.nii' ...
};

% Parameters for band-pass filtering.
num_vols_to_discard = 15; % Discard first 15 volumes
TR = 1; % Repetition Time in seconds. **ADJUST IF YOUR TR IS DIFFERENT**
low_cutoff_freq = 0.01; % Low cutoff for band-pass filter (Hz)
high_cutoff_freq = 0.08; % High cutoff for band-pass filter (Hz)
add_mean_back = 'Yes'; % Add the mean back to the time series after filtering.

% Parameters for memory management during filtering.
% Increase if you have more RAM, decrease if you run into memory errors.
CUTNUMBER = 10; 

% --- PROCESSING LOOP ---
% NOTE: The original code loops i=3:33. Ensure this range covers all your desired subjects.
for i = 3:Length_Names
    subject_id = FileNames{i};
    fprintf('--- Starting temporal preprocessing for subject: %s ---\n', subject_id);
    
    % Create a subject-specific output directory.
    subject_out_dir = fullfile(path_processed_fmri, subject_id);
    if ~exist(subject_out_dir, 'dir'), mkdir(subject_out_dir); end
    
    % Initialize matrices to store concatenated time series from all 4 runs.
    TS_all_runs = [];     % Whole brain
    ROI_L_all_runs = [];  % Left ROI
    ROI_R_all_runs = [];  % Right ROI

    for j = 1:numel(run_templates)
        % Load the smoothed CIFTI data.
        cifti_path = fullfile(path_smoothed_fmri, subject_id, run_templates{j});
        if ~exist(cifti_path, 'file')
            fprintf('WARNING: Smoothed file not found, skipping run %d for subject %s\n', j, subject_id);
            continue;
        end
        
        smooth_data = cifti_read(cifti_path);
        
        % Discard initial volumes.
        AllVolume = smooth_data.cdata(:, (num_vols_to_discard+1):end);
        [nDimVertex, nDimTimePoints] = size(AllVolume);
        
        % --- Band-Pass Filtering ---
        % The time series is transposed (time x vertices) for filtering.
        AllVolume = AllVolume';
        
        % Remove the mean before filtering.
        theMean = mean(AllVolume);
        AllVolume = AllVolume - repmat(theMean, [nDimTimePoints, 1]);
        
        % Process in segments to conserve memory.
        SegmentLength = ceil(size(AllVolume, 2) / CUTNUMBER);
        for iCut = 1:CUTNUMBER
            if iCut ~= CUTNUMBER
                Segment = (iCut - 1) * SegmentLength + 1 : iCut * SegmentLength;
            else
                Segment = (iCut - 1) * SegmentLength + 1 : size(AllVolume, 2);
            end
            AllVolume(:, Segment) = y_IdealFilter(AllVolume(:, Segment), TR, [low_cutoff_freq, high_cutoff_freq]);
        end
        
        % Add the mean back if specified.
        if strcmpi(add_mean_back, 'Yes')
            AllVolume = AllVolume + repmat(theMean, [nDimTimePoints, 1]);
        end
        
        % Transpose back to (vertices x time).
        ts_filtered = AllVolume';
        
        % --- ROI Extraction ---
        % !! IMPORTANT !!
        % The indices below are specific to a particular brain atlas and parcellation.
        % You MUST replace them with the correct vertex indices for your ROIs.
        ROI_L_indices = 160307:160889; % Example: Left Globus Pallidus
        ROI_R_indices = 160890:161388; % Example: Right Globus Pallidus
        
        ROI_L_ts = ts_filtered(ROI_L_indices, :);
        ROI_R_ts = ts_filtered(ROI_R_indices, :);
        
        % --- Normalization (z-scoring) and Concatenation ---
        % Transpose to (time x vertices) for normalization.
        ts_run = ts_filtered';
        ROI_L_ts_run = ROI_L_ts';
        ROI_R_ts_run = ROI_R_ts';
        
        % Detrend (remove mean) and scale to standard deviation of 1.
        ts_run = detrend(ts_run, 'constant'); ts_run = ts_run ./ repmat(std(ts_run), size(ts_run, 1), 1);
        ROI_L_ts_run = detrend(ROI_L_ts_run, 'constant'); ROI_L_ts_run = ROI_L_ts_run ./ repmat(std(ROI_L_ts_run), size(ROI_L_ts_run, 1), 1);
        ROI_R_ts_run = detrend(ROI_R_ts_run, 'constant'); ROI_R_ts_run = ROI_R_ts_run ./ repmat(std(ROI_R_ts_run), size(ROI_R_ts_run, 1), 1);
        
        % Concatenate this run's data with data from previous runs.
        TS_all_runs = [TS_all_runs; ts_run];
        ROI_L_all_runs = [ROI_L_all_runs; ROI_L_ts_run];
        ROI_R_all_runs = [ROI_R_all_runs; ROI_R_ts_run];
        
        fprintf('  - Finished run %d for subject %s\n', j, subject_id);
    end
    
    % Final normalization of the fully concatenated time series.
    TS_all_runs = detrend(TS_all_runs, 'constant'); TS_all_runs = TS_all_runs ./ repmat(std(TS_all_runs), size(TS_all_runs, 1), 1);
    ROI_L_all_runs = detrend(ROI_L_all_runs, 'constant'); ROI_L_all_runs = ROI_L_all_runs ./ repmat(std(ROI_L_all_runs), size(ROI_L_all_runs, 1), 1);
    ROI_R_all_runs = detrend(ROI_R_all_runs, 'constant'); ROI_R_all_runs = ROI_R_all_runs ./ repmat(std(ROI_R_all_runs), size(ROI_R_all_runs, 1), 1);
    
    % Save the final concatenated and preprocessed time series for the subject.
    save(fullfile(subject_out_dir, 'TS.mat'), 'TS_all_runs', '-v7.3');
    save(fullfile(subject_out_dir, 'ROI_L_TS.mat'), 'ROI_L_all_runs', '-v7.3');
    save(fullfile(subject_out_dir, 'ROI_R_TS.mat'), 'ROI_R_all_runs', '-v7.3');
    
    fprintf('Finished temporal preprocessing for subject: %s\n', subject_id);
end


%% ========================================================================
%  PART 3: FC AND SIMILARITY MATRIX CALCULATION
% =========================================================================
% This section calculates functional connectivity for each subject.
% The approach is to:
%   1. Use PCA to find the dominant whole-brain activity patterns (components).
%   2. Correlate the ROI time series with these principal components.
%   3. Create a vertex-wise similarity matrix for each ROI based on these
%      connectivity profiles.
% -------------------------------------------------------------------------

for i = 3:Length_Names
    subject_id = FileNames{i};
    subject_dir = fullfile(path_processed_fmri, subject_id);
    fprintf('--- Starting FC calculation for subject: %s ---\n', subject_id);

    % Load the preprocessed time series data.
    load(fullfile(subject_dir, 'TS.mat')); % Loads TS_all_runs
    load(fullfile(subject_dir, 'ROI_L_TS.mat')); % Loads ROI_L_all_runs
    load(fullfile(subject_dir, 'ROI_R_TS.mat')); % Loads ROI_R_all_runs

    % Perform PCA on the whole-brain time series to get principal components.
    % The 'scores' matrix contains the time course of each component.
    [~, scores] = pca(TS_all_runs);
    
    % Calculate the correlation between each ROI vertex and the brain-wide components.
    % This results in a (num_roi_vertices x num_components) connectivity matrix.
    % atanh() is the Fisher's z-transform, used to normalize correlation coefficients.
    FC_L = atanh(corr(ROI_L_all_runs, scores));
    FC_R = atanh(corr(ROI_R_all_runs, scores));
    
    % Save the connectivity matrices.
    save(fullfile(subject_dir, 'FC_L.mat'), 'FC_L', '-v7.3');
    save(fullfile(subject_dir, 'FC_R.mat'), 'FC_R', '-v7.3');
    
    % Calculate a vertex-wise similarity matrix from the FC profiles.
    % This measures how similarly any two vertices within the ROI are connected
    % to the rest of the brain.
    % NOTE: 'eta_squared' is a custom function. It likely computes the squared
    % correlation (R^2) between the connectivity profiles of each pair of vertices.
    FCS_L = eta_squared(FC_L); % Results in a (num_roi_vertices x num_roi_vertices) matrix
    FCS_R = eta_squared(FC_R);
    
    % Save the final similarity matrices for the subject.
    save(fullfile(subject_dir, 'FCS_L.mat'), 'FCS_L', '-v7.3');
    save(fullfile(subject_dir, 'FCS_R.mat'), 'FCS_R', '-v7.3');
    
    fprintf('Finished FC calculation for subject: %s\n', subject_id);
end


%% ========================================================================
%  PART 4: GROUP-LEVEL AGGREGATION
% =========================================================================
% This section loads the individual subject similarity matrices and computes
% a group average. It then combines this result with pre-computed DTI data.
% -------------------------------------------------------------------------

% --- SETUP ---
num_subjects = Length_Names - 2;
% !! NOTE !!: The dimensions for the pre-allocated matrices must match your ROI sizes.
% Please verify these numbers based on your atlas.
num_vertices_L = 583; % Example from original code, change to your data
num_vertices_R = 499; % Example from original code, change to your data

% Pre-allocate matrices to hold data from all subjects.
FC_similarity_L_all = zeros(num_vertices_L, num_vertices_L, num_subjects);
FC_similarity_R_all = zeros(num_vertices_R, num_vertices_R, num_subjects);

% --- AGGREGATION LOOP ---
for i = 3:Length_Names
    subject_id = FileNames{i};
    subject_dir = fullfile(path_processed_fmri, subject_id);
    
    fprintf('Aggregating data for subject: %s\n', subject_id);
    
    % Load the similarity matrices for the current subject.
    load(fullfile(subject_dir, 'FCS_L.mat')); % Loads FCS_L
    load(fullfile(subject_dir, 'FCS_R.mat')); % Loads FCS_R
    
    % Add the subject's data to the group matrix.
    FC_similarity_L_all(:, :, i-2) = FCS_L;
    FC_similarity_R_all(:, :, i-2) = FCS_R;
end

% Compute the average similarity matrix across all subjects.
FC_similarity_L_mean = mean(FC_similarity_L_all, 3);
FC_similarity_R_mean = mean(FC_similarity_R_all, 3);

% --- FINAL STEP: Combine with DTI data ---
% This step assumes you have pre-computed DTI-based similarity matrices
% loaded into the workspace (e.g., 'dti_similarity_L' and 'dti_similarity_R').

% Load your DTI data here if it's not already in the workspace.
% load('/path/to/your/dti_similarity_matrices.mat');

% Concatenate the DTI and fMRI similarity matrices for further analysis.
FC_DTI_combined_L = [dti_similarity_L_mean, FC_similarity_L_mean];
FC_DTI_combined_R = [dti_similarity_R_mean, FC_similarity_R_mean];

fprintf('--- Pipeline complete. Group average matrices are computed and combined with DTI data. ---\n');


%% Shared Nearest Neighbor Density Peak Clustering (SNNDPC)
%
% This script implements the SNNDPC algorithm, a variant of the Density Peak 
% Clustering (DPC) method. It is designed to identify cluster centers based on
% both high local density (rho) and large separation distance (delta). The "SNN"
% aspect incorporates a shared nearest neighbor similarity metric to better 
% handle clusters of varying shapes and densities.
%
% Pipeline:
%   1. Prepare data and set parameters.
%   2. Calculate core DPC metrics (rho, delta) using an SNN-based similarity.
%   3. Automatically identify cluster centers based on the gamma metric.
%   4. Assign all non-center points to clusters in a two-phase process.
%   5. (Optional) Visualize the results in a data-specific context (e.g., as a brain map).
%
% Dependencies:
%   - MATLAB (Statistics and Machine Learning Toolbox for 'pdist')
%   - (Optional) Tools for reading/writing specific data formats (e.g., NIfTI tools).
%
%
% Reference: This implementation is based on the concepts described in the
% original DPC paper by Rodriguez & Laio (2014) and incorporates SNN principles.
%

%% ========================================================================
%  1. DATA PREPARATION & PARAMETER SETUP
% =========================================================================

% --- USER-DEFINED PARAMETERS ---
% Load your data here. 'data' should be a matrix where each row is a data point
% and each column is a feature.
s = FC_DTI_combined_L; % Each hemisphere calculates separately.
data = s; % Example: s is a pre-loaded variable from your workspace

% NC: The number of clusters to identify.
num_clusters = your_cluster_count; % Example: 3, determine your cluster count based on decision graphs and other factors.

% Visualization flag and parameters (only used in Part 5).
enable_visualization = true; % Set to 'false' to skip the visualization part.
path_to_mask_file = '/path/to/your/mask/file.nii'; % Example for neuroimaging data
path_to_output_dir = '/path/to/your/output/directory';
output_filename = 'Clustering_Result.nii';
K = 30;% Nearest neighbors;
% --- END OF PARAMETERS ---

% Get the total number of data points.
N = size(data, 1);

% Optional: Normalize data features to a [0, 1] range.
% This can be important if features have vastly different scales.
% Uncomment the lines below if normalization is needed.
% data = (data - min(data)) ./ (max(data) - min(data));
% data(isnan(data)) = 0;


%% ========================================================================
%  2. CORE SNNDPC ALGORITHM: CALCULATE RHO AND DELTA
% =========================================================================
% This section calculates the two fundamental quantities for DPC:
%   - rho (ρ): Local density of each point.
%   - delta (δ): Minimum distance to a point with higher density.
% -------------------------------------------------------------------------

fprintf('Calculating pairwise distances...\n');
% --- Step 2a: Calculate initial pairwise distances ---
% 'cityblock' (Manhattan distance) is used here, but 'euclidean' is also common.
dist1 = squareform(pdist(data, 'cityblock'));

% --- Step 2b: Calculate SNN-based similarity (dist2) ---
% This is the key modification from standard DPC. Instead of using raw
% distance for density, we compute a similarity score based on shared neighbors.
fprintf('Calculating SNN similarity...\n');
[dist1Sort, dist1Order] = sort(dist1, 2);

% Pre-sort neighbor lists for efficient comparison.
dist1OrderSort = sort(dist1Order(:, 1:K), 2); 

dist2 = zeros(N);
sharedCount = zeros(N);
for p = 1:N
    % A clever optimization: check which neighbors of point p are also potential
    % neighbors of all other points further down the list.
    isNeighbor = dist1(p+1:N, dist1OrderSort(p, :)) <= dist1Sort(p+1:N, K);
    
    for o = p+1:N
        % Find shared neighbors between point p and o very efficiently.
        shared_neighbors = dist1OrderSort(p, isNeighbor(o-p, :));
        sharedCount(p, o) = length(shared_neighbors);
        
        % Calculate similarity only if p and o are in each other's K-neighborhood.
        % This ensures the connection is significant.
        if dist1(p, o) < min(dist1Sort(p, K+1), dist1Sort(o, K+1))
            % The similarity metric is based on the number of shared neighbors
            % normalized by the sum of their distances to p and o.
            dist2(p, o) = sharedCount(p, o)^2 / sum(dist1(p, shared_neighbors) + dist1(o, shared_neighbors));
        end
    end
end
% Make matrices symmetric.
dist2 = dist2 + dist2';
sharedCount = sharedCount + sharedCount';

% --- Step 2c: Calculate local density (rho) ---
% Rho is defined as the sum of SNN similarities to the K nearest neighbors.
fprintf('Calculating local density (rho)...\n');
dist2Sort = sort(dist2, 2, 'descend');
rho = sum(dist2Sort(:, 1:K), 2)';

% --- Step 2d: Calculate minimum distance to higher density point (delta) ---
fprintf('Calculating separation distance (delta)...\n');
delta = inf(1, N);
deltaSelect = zeros(1, N); % Stores the index of the point that defines delta
[~, rhoOrder] = sort(rho, 'descend');

% Pre-calculate sum of distances for a weighted delta calculation.
dist1SortSum = sum(dist1Sort(:, 1:K), 2);

for p = 2:N
    current_point_idx = rhoOrder(p);
    higher_density_indices = rhoOrder(1:p-1);
    
    % Calculate a weighted distance to all points with higher density.
    delta_candidates = dist1(current_point_idx, higher_density_indices) .* ...
                      (dist1SortSum(current_point_idx) + dist1SortSum(higher_density_indices)');
    
    % Find the minimum among these candidates.
    [delta(current_point_idx), min_idx] = min(delta_candidates);
    deltaSelect(current_point_idx) = higher_density_indices(min_idx);
end

% The point with the global maximum density has no higher-density point;
% its delta is defined as the maximum of all other delta values.
delta(rhoOrder(1)) = max(delta(delta ~= inf));

% --- Step 2e: Calculate the decision metric (gamma) ---
% Gamma is the product of rho and delta, combining both metrics.
% High gamma values are excellent candidates for cluster centers.
gamma = rho .* delta;
[sort_gamma, ~] = sort(gamma, 'descend');

% --- Step 2f: Generate decision plots ---
fprintf('Generating decision plots...\n');
figure('Name', 'DPC Decision Graph');
scatter(rho, delta, 30, 'filled', 'MarkerFaceAlpha', 0.7);
xlabel('\rho (Local Density)');
ylabel('\delta (Separation Distance)');
title('Decision Graph');

figure('Name', 'Gamma Plot');
scatter(1:N, sort_gamma, 30, 'filled');
xlabel('Point Index (Sorted by \gamma)');
ylabel('\gamma Value');
title('Gamma = \rho \times \delta');


%% ========================================================================
%  3. IDENTIFY CLUSTER CENTERS
% =========================================================================
% This section automatically selects cluster centers by finding the points with
% the highest gamma values.
% -------------------------------------------------------------------------
fprintf('Identifying cluster centers...\n');
% Initialize cluster assignments: -1 indicates an unassigned point.
cluster = -ones(1, N);

% Select the top 'num_clusters' points with the highest gamma values as centers.
center_indices = find(gamma >= sort_gamma(num_clusters));

% Assign a unique cluster ID (1, 2, 3...) to each center.
NC = length(center_indices);
cluster(center_indices) = 1:NC;


%% ========================================================================
%  4. ASSIGN NON-CENTER POINTS TO CLUSTERS
% =========================================================================
% This section assigns the remaining points using a two-phase approach.
% -------------------------------------------------------------------------

% --- Phase 4a: Assign high-confidence points via neighbor propagation ---
% This phase assigns points that are strongly connected to the cluster cores.
% It uses a queue-based approach (Breadth-First Search) starting from the centers.
fprintf('Assigning high-confidence points...\n');
queue = java.util.LinkedList(); % Using Java's LinkedList as a queue
for p = center_indices
    queue.add(p);
end

while ~queue.isEmpty()
    this_point = queue.remove();
    % Check the K-nearest neighbors of the current point.
    for next_point = dist1Order(this_point, 2:K+1)
        % Assign if the neighbor is unassigned AND shares many neighbors with the current point.
        if cluster(next_point) < 0 && sharedCount(this_point, next_point) >= K / 2
            cluster(next_point) = cluster(this_point);
            queue.add(next_point); % Add the newly assigned point to the queue to propagate its cluster
        end
    end
end

% --- Phase 4b: Assign remaining "halo" points ---
% This phase assigns any points left over from the first phase. It uses a
% voting system based on the clusters of each point's neighbors.
fprintf('Assigning remaining halo points...\n');
unassigned_indices = find(cluster < 0);
K_dynamic = K; % Use a dynamic K that can expand if needed

while ~isempty(unassigned_indices)
    % Create a "recognition matrix" for voting.
    recog = zeros(length(unassigned_indices), NC);
    for p = 1:length(unassigned_indices)
        current_point = unassigned_indices(p);
        % Check the neighbors of the unassigned point.
        for o = 2:K_dynamic + 1
            neighbor_cluster = cluster(dist1Order(current_point, o));
            % If the neighbor belongs to a cluster, cast a vote.
            if neighbor_cluster > 0
                recog(p, neighbor_cluster) = recog(p, neighbor_cluster) + 1;
            end
        end
    end
    
    % Find the highest vote count.
    max_votes = max(recog(:));
    
    if max_votes > 0
        % Find which points and clusters received the maximum votes.
        [points_to_assign_idx, clusters_to_assign] = find(recog == max_votes);
        % Assign the points to their new clusters.
        cluster(unassigned_indices(points_to_assign_idx)) = clusters_to_assign;
        % Update the list of unassigned points.
        unassigned_indices = find(cluster < 0);
    else
        % If no votes were cast (e.g., a point is too isolated), expand the
        % search neighborhood and try again.
        K_dynamic = K_dynamic + 1;
        fprintf('Expanding neighborhood search to K = %d\n', K_dynamic);
    end
end

% Ensure final output arrays have a consistent column vector shape.
if ~iscolumn(center_indices)
    center_indices = center_indices';
end
if ~iscolumn(cluster)
    cluster = cluster';
end
fprintf('Clustering complete. Found %d clusters.\n', NC);

%% Visualization
roiFile='E:\research\GP-sub\7T_result\GP_mask_7T_left.nii'; % Replace with the local address of your mask
%subcortex mask
[~,ins_roi]=read(roiFile); ind_roi=find(ins_roi);
GP_lh=zeros(113,136,113);
GP_lh(ind_roi)=cluster;
mat2nii(GP_lh,'E:\research\GP-sub\7T_result\GP_left3T.nii',size(GP_lh),32,'E:\research\GP-sub\7T_result\GP_mask_7T_left.nii');