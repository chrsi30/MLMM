function [xbin,ybinquantile] = binnedquantilesSB(x,y,quantile_value,binningInfo,FLAGlogX)

if nargin <= 3,
    binningInfo = 15;  % number bins
    FLAGlogX = 0;
elseif nargin<=4,
    FLAGlogX = 0;
end    

if length(binningInfo) == 1,
    % Bin equidistantly
    numbins = binningInfo;
    if FLAGlogX,
        bins = logspace(log10(min(x)), log10(max(x)), numbins);
    else
        bins = linspace(min(x), max(x), numbins);
    end
    
    try
        [n,bin] = histc(x, bins); %#ok<ASGLU>
    catch
        'error'
    end
    
    mu = NaN*zeros(size(bins));
    for k = [1:numbins], %#ok<NBRAK>
        ind = find(bin==k);
        if (~isempty(ind))
            mu(k) = quantile(y(ind),quantile_value);
        end
    end
    
    % Remove NaNs
    Z = [bins(:) mu(:)];
    Z(isnan(Z(:,2)),:) = [];
    
    % Assign outputs
    xbin = Z(:,1);
    ybinquantile = Z(:,2);
else
    % bin my mean binning value and look around range
    bins_mean = binningInfo{1};
    bins_lookaround = binningInfo{2};
    
    ybinquantile = [];
    for k=1:length(bins_mean),
        ix = find(x>=bins_mean(k)-bins_lookaround(k) & x<bins_mean(k)+bins_lookaround(k));
        ybinquantile(k) = quantile(y(ix),quantile_value);
    end
    xbin = bins_mean;
    ix = find(isnan(ybinquantile));
    xbin(ix) = [];
    ybinquantile(ix) = [];
end

