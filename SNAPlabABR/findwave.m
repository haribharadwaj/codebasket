function [peakval, peakind, troughval, troughind] = findwave(ref, y, whichwave)

switch whichwave
    % Here, we have to manually specify and peak and trough index for the
    % reference waveform
    case 0 % SP
        % Note trough occurs before peak for this one
        inds = [79, 70]; 
    case 1 % Wave I
        inds = [90, 100];
    case 5 % Wave 5
        inds = [156, 168];
    otherwise
        error('Unknown wave');
end


[~, iref, iy] = dtw(ref, y);

peakinds = iy(iref == inds(1));
troughinds = iy(iref == inds(2));
[peakval, ind1] = max(y(peakinds));
[troughval, ind2] = min(y(troughinds));

peakind = peakinds(ind1);
troughind = troughinds(ind2);