function [damage, power, level, depth] = getResults(fname)

load(fname);
c = 0;
cols = 'rgbyck';
for damage = 0:0.2:1
    w = [(1-damage)*0.15, (1-damage)*0.25, 0.6];
    mixedSR = w*pop_resp;
    c = c+1;
    [pxx,f] = pmtm(mixedSR,1.25,[],1/mean(diff(t_gated)));
    [f0,f0_ind] = min(abs(f - 100)); %#ok<ASGLU>
    P(c) = pxx(f0_ind);
    level = stimdb;
    depth = m*100;
end

damage = 100*(0:0.2:1);
power = P; %/P(1);
