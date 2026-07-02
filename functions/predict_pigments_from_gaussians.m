function [pigTbl, modelInfo] = predict_pigments_from_gaussians(G, model_file)
%PREDICT_PIGMENTS_FROM_GAUSSIANS Apply log-log MLR pigment algorithm.
%
% Model workbook sheets used:
%   TChla_coeffs              : TChla from agaus_676_newpeaks
%   all13_global4_coeffs      : global 4-Gaussian coefficients
%   HexBut_regime_coeffs      : low/high coefficients and TChla threshold
%
% Formula:
%   log10(pigment) = B0 + B1*log10(G434) + B2*log10(G496) +
%                    B3*log10(G524) + B4*log10(G676)

if ~isfile(model_file)
    error('Model coefficient workbook not found: %s', model_file);
end

requiredG = {'agaus_434_newpeaks','agaus_496_newpeaks','agaus_524_newpeaks','agaus_676_newpeaks'};
for i = 1:numel(requiredG)
    if ~ismember(requiredG{i}, G.Properties.VariableNames)
        error('Missing required Gaussian predictor: %s', requiredG{i});
    end
end

n = height(G);
Xraw = [G.agaus_434_newpeaks, G.agaus_496_newpeaks, G.agaus_524_newpeaks, G.agaus_676_newpeaks];
Xraw(~isfinite(Xraw) | Xraw <= 0) = NaN;
X = log10(Xraw);
Xdesign = [ones(n,1), X];
validX = all(isfinite(X), 2);

%% TChla
TchlCoeff = readtable(model_file, 'Sheet','TChla_coeffs', 'VariableNamingRule','preserve');
Bchl = [TchlCoeff.Intercept(1); TchlCoeff.Slope(1)];
G676 = G.agaus_676_newpeaks;
TChla = nan(n,1);
m676 = isfinite(G676) & G676 > 0;
TChla(m676) = 10.^(Bchl(1) + Bchl(2).*log10(G676(m676)));

pigTbl = table(TChla, 'VariableNames', {'TChla'});

%% Global coefficients for all pigments
Cglob = readtable(model_file, 'Sheet','all13_global4_coeffs', 'VariableNamingRule','preserve');
for i = 1:height(Cglob)
    pig = string(Cglob.Pigment(i));
    B = [Cglob.Intercept(i); Cglob.G434(i); Cglob.G496(i); Cglob.G524(i); Cglob.G676(i)];

    y = nan(n,1);
    y(validX) = 10.^(Xdesign(validX,:) * B);
    pigTbl.(pig) = y;
end

%% Replace Hex and But with two-regime predictions
Creg = readtable(model_file, 'Sheet','HexBut_regime_coeffs', 'VariableNamingRule','preserve');
TChla_threshold = Creg.TChla_threshold(1);

for i = 1:height(Creg)
    pig = string(Creg.Pigment(i));

    Blow = [Creg.low_Intercept(i); Creg.low_G434(i); Creg.low_G496(i); Creg.low_G524(i); Creg.low_G676(i)];
    Bhigh = [Creg.high_Intercept(i); Creg.high_G434(i); Creg.high_G496(i); Creg.high_G524(i); Creg.high_G676(i)];

    y = nan(n,1);
    mLow = validX & isfinite(TChla) & TChla > 0 & TChla <= TChla_threshold;
    mHigh = validX & isfinite(TChla) & TChla > TChla_threshold;

    y(mLow) = 10.^(Xdesign(mLow,:) * Blow);
    y(mHigh) = 10.^(Xdesign(mHigh,:) * Bhigh);

    pigTbl.(pig) = y;
end

modelInfo = struct();
modelInfo.model_file = model_file;
modelInfo.predictors = requiredG;
modelInfo.TChla_threshold_for_HexBut = TChla_threshold;
modelInfo.formula = 'log10(y) = intercept + coefficients*log10(agaus)';
modelInfo.n_valid_predictor_rows = nnz(validX);
end
