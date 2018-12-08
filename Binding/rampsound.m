function y = rampsound(x,fs,risetime)
% Function to ramp a sound file using a dpss ramp
% USAGE:
% y = rampsound(x,fs,risetime)
%
% risetime in seconds, fs in Hz
% Hari Bharadwaj

Nramp = ceil(fs*risetime*2)+1;
w = dpss(Nramp,1,1);

w = w - w(1);
w = w/max(w);
sz = size(x);
if( sz(1) < sz(2) )
    nx = sz(2);
    vert = 1;
else
    nx = sz(1);
    vert = 0;
end

half = ceil(Nramp/2);
wbig = [w(1:half); ones(nx - 2*half,1); w((end-half+1):end)];



if(nx == numel(x))
    if(vert)
        y = x.*wbig';
    else
        y = x.*wbig;
    end
else
    if(vert)
        y = x.*repmat(wbig, 1, 2)';
    else
        y = x.*repmat(wbig, 1, 2);
    end
end

