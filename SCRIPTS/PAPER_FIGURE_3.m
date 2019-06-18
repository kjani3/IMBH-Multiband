clear , clc


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

Mdet=logspace(1,7, 100); % from 10^1 to 10^7 solar mass
Ddet=100; % Mpc






%% (3a) Epoch-4: LISA

Detector = 'LISA';
timeline =4; 
rho_threshold  = 8;

Msrc_LISA = zeros(numel(Mdet),1);
dAvg_LISA = zeros(numel(Mdet),1);
zAvg_LISA = zeros(numel(Mdet),1);

for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    [fmin, fmax] = LISA_CBC_FMIN( q, Mdet(i), timeline);
    [rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dAvg_LISA(i), zAvg_LISA(i)] = AVG_DIST(rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list);

    Msrc_LISA(i) = Mdet(i)./(1+zAvg_LISA(i));    

end



%% (3b) Epoch-3: Einstein Telescope

Detector = 'ET';
rho_threshold  = 8;

Msrc_Ground = zeros(numel(Mdet),1);
dAvg_Ground = zeros(numel(Mdet),1);
zAvg_Ground = zeros(numel(Mdet),1);


for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    fmin = 1; fmax = 5e3;
    [rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dAvg_Ground(i), zAvg_Ground(i)] = AVG_DIST(rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list);

    Msrc_Ground(i) = Mdet(i)./(1+zAvg_Ground(i));    

end
    
% Multiband Distance

MMulti = logspace(1, log10(2e5), 100); % from 10^1 to 10^10 solar mass

[dMulti] = MULTIBAND_DISTANCE(MMulti, dAvg_LISA, Msrc_LISA, dAvg_Ground, Msrc_Ground);


set(gca,'fontsize',14)
area(MMulti, dMulti'./1e3,'LineStyle','-','LineWidth',1); hold on
set (gca, 'Yscale', 'log');
set (gca, 'Xscale', 'log');
xlim([10, 3e4])
ylim([1e-2, 1e2])
alpha(0.7)
grid on
axis square

%% (3b) Epoch-3: Cosmic Explorer

Detector = 'CE';
rho_threshold  = 8;

Msrc_Ground = zeros(numel(Mdet),1);
dAvg_Ground = zeros(numel(Mdet),1);
zAvg_Ground = zeros(numel(Mdet),1);


for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    fmin = 5; fmax = 2e3;
    [rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dAvg_Ground(i), zAvg_Ground(i)] = AVG_DIST(rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list);

    Msrc_Ground(i) = Mdet(i)./(1+zAvg_Ground(i));    

end
    
% Multiband Distance

MMulti = logspace(1, log10(2e5), 100); % from 10^1 to 10^10 solar mass

[dMulti] = MULTIBAND_DISTANCE(MMulti, dAvg_LISA, Msrc_LISA, dAvg_Ground, Msrc_Ground);
area(MMulti, dMulti'./1e3,'LineStyle','-','LineWidth',1); hold on

%% (3b) Epoch-3: Voyager

Detector = 'VOY';
rho_threshold  = 8;

Msrc_Ground = zeros(numel(Mdet),1);
dAvg_Ground = zeros(numel(Mdet),1);
zAvg_Ground = zeros(numel(Mdet),1);


for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    fmin = 5; fmax = 2e3;
    [rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dAvg_Ground(i), zAvg_Ground(i)] = AVG_DIST(rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list);
    Msrc_Ground(i) = Mdet(i)./(1+zAvg_Ground(i));    

end
    
% Multiband Distance

MMulti = logspace(1, log10(2e5), 100); % from 10^1 to 10^10 solar mass

[dMulti] = MULTIBAND_DISTANCE(MMulti, dAvg_LISA, Msrc_LISA, dAvg_Ground, Msrc_Ground);


area(MMulti, dMulti'./1e3,'LineStyle','-','LineWidth',1); hold on


%% (5) Plot


legend('Network: LISA + ET', 'Network: LISA + CE', 'Network: LISA + Voy')
xlabel('Total Mass in Source Frame [$M_\odot$]','Interpreter','latex','FontName','TimesNewRoman','FontSize',14)
ylabel('Multiband Detection Radius [Gpc]','Interpreter','latex','FontName','TimesNewRoman','FontSize',14)
alpha(0.7)



