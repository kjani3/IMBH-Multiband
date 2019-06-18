function [ galactic_noise ] = LISA_GAL_NOISE( freq, timeline)
% LISA_GAL_NOISE computes galactic background noise in LISA. 
% Reference: arXiv 1803.01944

if timeline == 0.5
    alpha = 0.133;
    beta = 243;
    kappa = 482;
    gamma = 917;
    fk = 0.00258;

elseif timeline == 1
    alpha = 0.171;
    beta = 299;
    kappa = 611;
    gamma = 1340;
    fk = 0.00215;
    
elseif timeline == 2
    alpha = 0.165;
    beta = 299;
    kappa = 611;
    gamma = 1340;
    fk = 0.00173;
    
elseif timeline == 4
    alpha = 0.138;
    beta = -221;
    kappa = 521;
    gamma = 1680;
    fk = 0.00113;
end

 A = 9e-45;
 
 galactic_noise  = A.*freq.^(-7/3).*exp(-freq.*alpha + beta.*freq.*sin(kappa.*freq)).*(1 + tanh(gamma.*(fk - freq)));


end

