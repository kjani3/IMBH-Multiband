function [dMulti] = MULTIBAND_DISTANCE(MMulti, dAvg_LISA, Msrc_LISA, dAvg_Ground, Msrc_Ground) 
% MULTIBAND_DISTANCE computes distance for the multiband network between LISA and ground detector.
% The input are MMulti (range of masses in source frame to compute
% distahce), angle-average and source masses in LISA and ground detector


dMulti = zeros(numel(MMulti),1);

d_LISA = interp1(Msrc_LISA, dAvg_LISA, MMulti);
d_Ground = interp1(Msrc_Ground, dAvg_Ground, MMulti);

j=0;
for i = 1:numel(MMulti)
    j = j+1;
    if d_LISA(i) >= d_Ground(i)
        temp_dist(j) = d_Ground(i);        
    else
        temp_dist(j) = d_LISA(i);
    end
    
    dMulti(j) = temp_dist(j);
    
end

end

