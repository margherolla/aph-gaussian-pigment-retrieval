function [T, aph, wl, aph_sd, wl_sd, metaTbl] = load_aph_excel(input_file)
%LOAD_APH_EXCEL Read station aph(lambda) spectra and uncertainties from Excel.
%
% Required aph column names:
%   aph_400, aph_401, ...
%   aph_(400), aph_(401), ...
%   aph_400p5 for decimal wavelengths such as 400.5
%
% Optional uncertainty column names:
%   aph_sd_400, aph_sd_401, ...
%   aph_sd_(400), aph_sd_(401), ...
%   aph_sd_400p5 for decimal wavelengths such as 400.5
%
% Output:
%   T       = full imported table
%   aph     = nStations x nWavelengths matrix
%   wl      = wavelength vector for aph
%   aph_sd  = nStations x nWavelengths uncertainty matrix
%             If aph_sd columns are missing, this is filled with NaN.
%   wl_sd   = wavelength vector for aph_sd, matched to wl
%   metaTbl = all non-aph and non-aph_sd columns

if ~isfile(input_file)
    error('Input file not found: %s', input_file);
end

T = readtable(input_file, 'VariableNamingRule','preserve');
varNames = string(T.Properties.VariableNames);

isAph   = false(size(varNames));
isAphSd = false(size(varNames));
wl_all  = nan(size(varNames));
wl_sd_all = nan(size(varNames));

for i = 1:numel(varNames)
    v = varNames(i);

    % Examples: aph_443, aph_(443), aph_443.5, aph_443p5
    tokAph = regexp(v, '^aph_\(?([0-9]+(?:[\._p][0-9]+)?)\)?$', 'tokens', 'once');

    % Examples: aph_sd_443, aph_sd_(443), aph_sd_443.5, aph_sd_443p5
    tokSd = regexp(v, '^aph_sd_\(?([0-9]+(?:[\._p][0-9]+)?)\)?$', 'tokens', 'once');

    if ~isempty(tokAph)
        wtxt = tokAph{1};
        wtxt = replace(wtxt, 'p', '.');
        wtxt = replace(wtxt, '_', '.');
        wnum = str2double(wtxt);
        if isfinite(wnum)
            isAph(i) = true;
            wl_all(i) = wnum;
        end
    elseif ~isempty(tokSd)
        wtxt = tokSd{1};
        wtxt = replace(wtxt, 'p', '.');
        wtxt = replace(wtxt, '_', '.');
        wnum = str2double(wtxt);
        if isfinite(wnum)
            isAphSd(i) = true;
            wl_sd_all(i) = wnum;
        end
    end
end

if ~any(isAph)
    error(['No aph wavelength columns found. Use names like aph_400 or aph_(400). ', ...
           'Current columns are: %s'], strjoin(varNames, ', '));
end

aphVars = varNames(isAph);
wl = wl_all(isAph);

% Sort wavelengths from blue to red.
[wl, ord] = sort(wl(:));
aphVars = aphVars(ord);

aph = T{:, aphVars};
aph = double(aph);
aph(~isfinite(aph)) = NaN;
aph(aph < 0) = NaN;

% Read aph uncertainty columns, if present, and align them to wl.
aph_sd = nan(size(aph));
wl_sd = wl;

if any(isAphSd)
    sdVars = varNames(isAphSd);
    wl_sd_found = wl_sd_all(isAphSd);

    [wl_sd_found, ordSd] = sort(wl_sd_found(:));
    sdVars = sdVars(ordSd);

    aph_sd_found = double(T{:, sdVars});
    aph_sd_found(~isfinite(aph_sd_found)) = NaN;
    aph_sd_found(aph_sd_found < 0) = NaN;

    % Match each aph wavelength to the corresponding aph_sd wavelength.
    for iw = 1:numel(wl)
        j = find(abs(wl_sd_found - wl(iw)) < 1e-9, 1);
        if ~isempty(j)
            aph_sd(:, iw) = aph_sd_found(:, j);
        end
    end

    missingSd = all(isnan(aph_sd), 1);
    if any(missingSd)
        warning('Missing aph_sd columns for %d aph wavelengths. Those uncertainties are set to 5%% of aph.', nnz(missingSd));
        aph_sd(:, missingSd) = 0.05 .* abs(aph(:, missingSd));
    end
else
    warning('No aph_sd wavelength columns found. Using 5%% relative uncertainty as fallback.');
    aph_sd = 0.05 .* abs(aph);
end

% Avoid exact zero uncertainty where aph exists, because decomposition routines
% often use uncertainty as a weight.
mBadSd = (~isfinite(aph_sd) | aph_sd <= 0) & isfinite(aph) & aph > 0;
aph_sd(mBadSd) = 0.05 .* abs(aph(mBadSd));

metaTbl = T(:, ~(isAph | isAphSd));
end
