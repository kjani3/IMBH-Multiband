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

Mdet=logspace(1,10, 100); % from 10^1 to 10^10 solar mass
Ddet=100; % Mpc


%% (3a) Epoch-4: LISA

Detector = 'LISA';
timeline =4; 
rho_threshold  = 8;

Msrc = zeros(numel(Mdet),1);
dHor = zeros(numel(Mdet),1);
zHor = zeros(numel(Mdet),1);

for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    [fmin, fmax] = LISA_CBC_FMIN( q, Mdet(i), timeline);
    [ rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dHor(i), zHor(i)] = HOR_DIST( rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list );
    
    Msrc(i) = Mdet(i)./(1+zHor(i));    

end
    

figure (1)
set(gca,'fontsize',14)
area(Msrc, zHor,'LineStyle','-','LineWidth',1); hold on
set (gca, 'Yscale', 'log');
set (gca, 'Xscale', 'log');
xlim([10, 1e9])
ylim([4e-2, 2e3])


%% (3b) Epoch-3: Cosmic Explorer

Detector = 'CE';
rho_threshold  = 8;

Msrc = zeros(numel(Mdet),1);
dHor = zeros(numel(Mdet),1);
zHor = zeros(numel(Mdet),1);


for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    fmin = 5; fmax = 2e3;
    [ rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dHor(i), zHor(i)] = HOR_DIST( rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list );
    Msrc(i) = Mdet(i)./(1+zHor(i));    

end
    

area(Msrc, zHor,'LineStyle','-','LineWidth',1); hold on

%% (3b) Epoch-3: Einstein Telescope

Detector = 'ET';
rho_threshold  = 8;

Msrc = zeros(numel(Mdet),1);
dHor = zeros(numel(Mdet),1);
zHor = zeros(numel(Mdet),1);


for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    fmin = 1; fmax = 2e3;
    [ rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dHor(i), zHor(i)] = HOR_DIST( rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list );

    Msrc(i) = Mdet(i)./(1+zHor(i));    

end
    
area(Msrc, zHor,'LineStyle','-','LineWidth',1); hold on



%% (3c) Epoch-2: Voyager

Detector = 'VOY';
rho_threshold  = 8;

Msrc = zeros(numel(Mdet),1);
dHor = zeros(numel(Mdet),1);
zHor = zeros(numel(Mdet),1);


for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    fmin = 5; fmax = 2e3;
    [ rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dHor(i), zHor(i)] = HOR_DIST( rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list );

    Msrc(i) = Mdet(i)./(1+zHor(i));    

end
    

area(Msrc, zHor,'LineStyle','-','LineWidth',1); hold on

 
%% (3c) Epoch-2: A+

Detector = 'A+';
rho_threshold  = 8;

Msrc = zeros(numel(Mdet),1);
dHor = zeros(numel(Mdet),1);
zHor = zeros(numel(Mdet),1);


for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    fmin = 10; fmax = 2e3;
    [ rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dHor(i), zHor(i)] = HOR_DIST( rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list );

    Msrc(i) = Mdet(i)./(1+zHor(i));    

end
    

area(Msrc, zHor,'LineStyle','-','LineWidth',1); hold on

%% (3d) Epoch-1: aLIGO

Detector = 'aLIGO';
rho_threshold  = 8;

Msrc = zeros(numel(Mdet),1);
dHor = zeros(numel(Mdet),1);
zHor = zeros(numel(Mdet),1);


for i = 1:numel(Mdet)
    
    [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
    fmin = 10; fmax = 2e3;
    [ rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
    [dHor(i), zHor(i)] = HOR_DIST( rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list );

    Msrc(i) = Mdet(i)./(1+zHor(i));    

end
    

area(Msrc, zHor,'LineStyle','-','LineWidth',1); hold on

%%


alpha(0.6)
grid on
legend('Epoch-4: LISA', 'Epoch-3: CE', 'Epoch-3: ET', 'Epoch-2: Voy','Epoch-2: A+', 'Epoch-1: aLIGO')
xlabel('Total Mass in Source Frame [$M_\odot$]','Interpreter','latex','FontName','TimesNewRoman','FontSize',20)
ylabel('Detection Radius [$z$]','Interpreter','latex','FontName','TimesNewRoman','FontSize',20)



