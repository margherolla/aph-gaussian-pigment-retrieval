function [clusterTbl, ratioTbl] = cluster_pigments_stations(pigTbl, nPigmentClusters, fig_dir)
%CLUSTER_PIGMENTS_STATIONS Cluster pigments using pigment/sum-pigments fractions.
%
% This function:
%   1) computes pigment fractions:
%          Rnorm = pigment / sum(all pigments)
%   2) keeps only rows where all pigment fractions are finite and positive
%   3) clusters pigments using:
%          distance = correlation
%          linkage  = ward

if nargin < 2 || isempty(nPigmentClusters), nPigmentClusters = 3; end
if nargin < 3 || isempty(fig_dir), fig_dir = 'figures'; end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

out_dir = fullfile(fig_dir, '..', 'outputs');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

% Remove TChla if present, because clustering is based only on pigment fractions
pigNames = setdiff(pigTbl.Properties.VariableNames, {'TChla'}, 'stable');
n = height(pigTbl);

%% ------------------------------------------------------------------------
% 1) Pigment / sum(pigments)
%% ------------------------------------------------------------------------
P = nan(n, numel(pigNames));

for i = 1:numel(pigNames)
    P(:,i) = pigTbl.(pigNames{i});
end

validRows = all(isfinite(P) & P > 0, 2);

sumPig = sum(P, 2, 'omitnan');

Rnorm = nan(size(P));
mNorm = validRows & sumPig > 0;

Rnorm(mNorm,:) = P(mNorm,:) ./ sumPig(mNorm);

fprintf('\nHCA rows kept (pigment/sum pigments): %d / %d\n', nnz(mNorm), n);

ratioTbl = array2table(Rnorm, ...
    'VariableNames', strcat('frac_', pigNames, '_to_sumPigments'));

writetable(ratioTbl, ...
    fullfile(out_dir, 'pigment_sumPigment_fractions.csv'));

clusterTbl = table();
clusterTbl.ValidForClustering = mNorm;
clusterTbl.PigmentCluster = nan(n,1);

if nnz(mNorm) < 5
    warning('Too few valid rows for pigment clustering. Skipping clustering.');
    return
end

%% ------------------------------------------------------------------------
% 2) Pigment clustering
%% ------------------------------------------------------------------------
Rp = Rnorm(mNorm,:);

Tp = Rp';

Dp = pdist(Tp,'correlation');
Zp = linkage(Dp,'ward');
[coph_p, ~] = cophenet(Zp, Dp);

fprintf('Pigment cophenetic correlation = %.3f\n', coph_p);

Tpig = cluster(Zp, 'maxclust', nPigmentClusters);

pigClusterTbl = table(pigNames(:), Tpig(:), ...
    'VariableNames', {'Pigment','PigmentCluster'});

writetable(pigClusterTbl, fullfile(out_dir, 'pigment_clusters.csv'));

%% ------------------------------------------------------------------------
% 3) Dendrogram
%% ------------------------------------------------------------------------
cluster_colors_pig = [
    0.00 0.75 0.75;
    1.00 0.50 0.00;
    0.50 0.50 0.50;
    0.45 0.45 1.00;
    0.60 0.30 0.80
];

label_colors = zeros(numel(pigNames), 3);

for ic = 1:nPigmentClusters
    idx = Tpig == ic;
    col = cluster_colors_pig(min(ic, size(cluster_colors_pig,1)), :);
    label_colors(idx,:) = repmat(col, sum(idx), 1);
end

fig1 = figure('Color','w', 'Position',[200 200 900 700]);

dendrogram(Zp, 'Labels', pigNames, 'ColorThreshold', 'default');

xlabel('Pigments', 'FontSize', 20, 'FontWeight','bold');
ylabel('Distance', 'FontSize', 20, 'FontWeight','bold');

hold on
if nPigmentClusters == 3
    yline(1.5, 'k--', 'LineWidth', 2);
end
hold off

ax = gca;

D = findall(ax, 'Type', 'Line');
for i = 1:length(D)
    D(i).Color = [0 0 0];
    D(i).LineWidth = 2;
end

xticklabels_now = ax.XTickLabel;

for idx = 1:length(xticklabels_now)
    pigment_name = xticklabels_now{idx};
    pigment_idx = find(strcmp(pigNames, pigment_name), 1);

    if ~isempty(pigment_idx)
        pc = label_colors(pigment_idx,:);
        ax.XTickLabel{idx} = ...
            ['\color[rgb]{' num2str(pc) '}' pigment_name];
    end
end

set(ax, 'XTickLabelRotation', 45, 'FontSize', 18);
pbaspect([2.5 2 2]);

exportgraphics(fig1, fullfile(fig_dir,'pigment_dendrogram.png'), ...
    'Resolution', 300);

%% ------------------------------------------------------------------------
% 4) Metadata
%% ------------------------------------------------------------------------
clusterTbl.PigmentClusteringMethod(:) = ...
    {'R = pigment/sum(pigments); pdist correlation; ward linkage'};

clusterTbl.PigmentCophenetic(:) = coph_p;

end
