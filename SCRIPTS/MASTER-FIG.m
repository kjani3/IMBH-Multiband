clear all , clc


%% (1) GENERATE REDSHIFT VS. LUMINOSITY DISTANCE

CosmoPars = [67.66; 0.3111; 0.6889]; 
% PLANCK 2018 Table 2 (TT, TE, EE +LowE + Lensing + BAO) [H OmegM OmegLamba]

z_list =logspace(-4, log10(5e3), 5000);
[DL_list, ~, ~]  = lum_dist(z_list, CosmoPars );
DL_list = DL_list./1e6; %Get in Mpc


%% (2) LOAD WAVEFORM

wf_dir = importdata('../WF/f1_q1_M50_D100_i0_IMRHM.dat');
M0 = 50; % wf total-mass in Msun
D0 = 100; % wf distance in Mpc
q=1;  % mass-ratio

wf_data = downsample(wf_dir.data, 1);
f0 = wf_data(2:end,1);
hpf0 = abs((wf_data(2:end,2) + sqrt(-1).*wf_data(2:end,3)));

wf_fmin = 1e-5;
[hpf_fit, f_fit] = FIT_INSPIRAL(wf_fmin, hpf0, f0  );



%% (3) SPECIFY DETECTOR & SNR-THRESHOLD

Detector = 'LISA';
timeline =4; 
%[Sn, f_Sn, ~, ~] = GET_NOISE(Detector);
rho_threshold  = 8;


%% (4) COMPUTE SNR FOR DETECTOR MASSES

Mdet=logspace(1,10, 500); % from 10^1 to 10^10 solar mass
Ddet=100; % Mpc

for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    [fmin, fmax] = LISA_CBC_FMIN( q, Mdet(i), timeline);
    [ rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dHor(i), zHor(i)] = HOR_DIST( rho_opt, rho_threshold, Ddet, z_list, DL_list );

    Msrc(i) = Mdet(i)./(1+zHor(i));    

end
    


%% (5) PLOT 

loglog(Msrc, zHor, 'LineWidth',2)



%%
