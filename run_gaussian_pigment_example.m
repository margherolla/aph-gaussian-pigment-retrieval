%% run_gaussian_pigment_example.m
% Example workflow:
% 1) read station aph(lambda) spectra from Excel
% 2) run Gaussian decomposition using spectral_decomp_aph.m
% 3) derive TChla and diagnostic pigments using MLR coefficients
% 4) cluster pigments and stations from pigment/TChla ratios
% 5) export tables and simple figures

clear; clc; close all;

%% ------------------------- USER SETTINGS -------------------------------
input_file = fullfile('data','example_aph_stations.xlsx');
model_file = fullfile('models','Pigment_MLR_final_14vars_HexBut2reg.xlsx');

% Folder containing spectral_decomp_aph.m.
% If you copy spectral_decomp_aph.m into /functions, this line is enough.
addpath('functions');

out_dir = 'outputs';
fig_dir = 'figures';
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

peak_list = [410 434 451 474 496 524 560 594 634 658 676 700];

%% ------------------------- LOAD aph EXCEL ------------------------------
% Expected aph columns: aph_400, aph_401, ... OR aph_(400), aph_(401), ...
% Expected uncertainty columns: aph_sd_400, ... OR aph_sd_(400), ...
% Optional metadata columns: Station, Latitude, Longitude, Date, Time, etc.
[dataTbl, aph, wl, aphunc, wlunc, metaTbl] = load_aph_excel(input_file);

fprintf('Loaded %d stations and %d aph wavelengths from %s\n', ...
    size(aph,1), numel(wl), input_file);

%% ------------------------- GAUSSIAN DECOMPOSITION ----------------------
if exist('spectral_decomp_aph','file') ~= 2
    error(['spectral_decomp_aph.m was not found. ', ...
           'Copy it into the functions folder or add its folder with addpath().']);
end

% aphunc and wlunc are read from the same Excel file using aph_sd_(wl) columns.
% If some uncertainty columns are missing, load_aph_excel fills them with
% a simple 5 percent relative uncertainty fallback.
acs = 0;

[amps_newpeaks, compspec_newpeaks, sumspec_newpeaks] = ...
    spectral_decomp_aph(wl, aph, wlunc, aphunc, acs);

if size(amps_newpeaks,2) < numel(peak_list)
    error('spectral_decomp_aph returned fewer Gaussian amplitudes than expected. Check peak_list.');
end

G = array2table(amps_newpeaks(:,1:numel(peak_list)), ...
    'VariableNames', compose('agaus_%d_newpeaks', peak_list));

%% ------------------------- PREDICT PIGMENTS ----------------------------
[pigTbl, modelInfo] = predict_pigments_from_gaussians(G, model_file);

%% ------------------------- CLUSTERING ----------------------------------
nPigmentClusters = 3;

[clusterTbl, ratioTbl] = cluster_pigments_stations( ...
    pigTbl, nPigmentClusters, fig_dir);

%% ------------------------- EXPORT RESULTS ------------------------------
resultsTbl = [metaTbl, G, pigTbl, clusterTbl, ratioTbl];

writetable(resultsTbl, fullfile(out_dir,'gaussian_pigment_predictions.xlsx'));
writetable(resultsTbl, fullfile(out_dir,'gaussian_pigment_predictions.csv'));

save(fullfile(out_dir,'gaussian_pigment_results.mat'), ...
    'dataTbl','metaTbl','aph','wl','G','pigTbl','ratioTbl','clusterTbl', ...
    'amps_newpeaks','compspec_newpeaks','sumspec_newpeaks','modelInfo');

fprintf('\nDone. Results saved in %s and figures saved in %s.\n', out_dir, fig_dir);