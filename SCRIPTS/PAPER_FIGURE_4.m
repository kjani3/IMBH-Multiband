clear , clc


%% (1) GENERATE REDSHIFT VS. LUMINOSITY DISTANCE

CosmoPars = [67.66; 0.3111; 0.6889]; 
% PLANCK 2018 Table 2 (TT, TE, EE +LowE + Lensing + BAO) [H OmegM OmegLamba]

z_list =logspace(-4, log10(5e3), 5000);
[DL_list, ~, ~]  = lum_dist(z_list, CosmoPars );
DL_list = DL_list./1e6; %Get in Mpc


%% (2) LOAD WAVEFORM OF VARRYING MASS-RATIO (qq)

wf_q1 = importdata('../WF/q-series/f5_q1_M50_D100_i0_IMRHM.dat');       qq(1) = 1;
wf_q1p5 = importdata('../WF/q-series/f5_q1.5_M50_D100_i0_IMRHM.dat');   qq(2) = 1.5;
wf_q2 = importdata('../WF/q-series/f5_q2_M50_D100_i0_IMRHM.dat');       qq(3) = 2;
wf_q3 = importdata('../WF/q-series/f5_q3_M50_D100_i0_IMRHM.dat');       qq(4) = 3;
wf_q4 = importdata('../WF/q-series/f5_q4_M50_D100_i0_IMRHM.dat');       qq(5) = 4;
wf_q5 = importdata('../WF/q-series/f5_q5_M50_D100_i0_IMRHM.dat');       qq(6) = 5;
wf_q6 = importdata('../WF/q-series/f5_q6_M50_D100_i0_IMRHM.dat');       qq(7) = 6;
wf_q7 = importdata('../WF/q-series/f5_q7_M50_D100_i0_IMRHM.dat');       qq(8) = 7;
wf_q8 = importdata('../WF/q-series/f5_q8_M50_D100_i0_IMRHM.dat');       qq(9) = 8;
wf_q9 = importdata('../WF/q-series/f5_q9_M50_D100_i0_IMRHM.dat');       qq(10) = 9;
wf_q10 = importdata('../WF/q-series/f5_q10_M50_D100_i0_IMRHM.dat');     qq(11) = 10;
wf_q15 = importdata('../WF/q-series/f5_q15_M50_D100_i0_IMRHM.dat');     qq(12) = 15;
wf_q20 = importdata('../WF/q-series/f5_q20_M50_D100_i0_IMRHM.dat');     qq(13) = 20;
wf_q30 = importdata('../WF/q-series/f5_q30_M50_D100_i0_IMRHM.dat');     qq(14) = 30;
wf_q40 = importdata('../WF/q-series/f5_q40_M50_D100_i0_IMRHM.dat');     qq(15) = 40;
wf_q50 = importdata('../WF/q-series/f5_q50_M50_D100_i0_IMRHM.dat');     qq(16) = 50;

wf_list = {wf_q1, wf_q1p5, wf_q2, wf_q3, wf_q4, wf_q5, wf_q6, wf_q7, wf_q8, wf_q9, wf_q10, wf_q15, wf_q20, wf_q30, wf_q40, wf_q50};

eta = qq./(1+qq).^2;
M0 = 50; % wf total-mass in Msun
D0 = 100; % wf distance in Mpc

%% (3) CREATE VARIABLES

Mdet=logspace(1,7, 500); % from 10^1 to 10^7 solar mass
Ddet=100; % Mpc
rho_threshold  = 8;

Msrc_LISA = zeros(numel(Mdet),numel(wf_list));
dAvg_LISA = zeros(numel(Mdet),numel(wf_list));
zAvg_LISA = zeros(numel(Mdet),numel(wf_list));

Msrc_Ground = zeros(numel(Mdet),numel(wf_list));
dAvg_Ground = zeros(numel(Mdet),numel(wf_list));
zAvg_Ground = zeros(numel(Mdet),numel(wf_list));

MMulti = logspace(log10(50), log10(2e4), 500); % from 50 to 2*10^5 solar mass
dMulti = zeros(numel(MMulti),numel(wf_list));


%% (4) COMPUTE MULTIBAND DISTANCE


for n = 1:numel(wf_list)

    wf_data = downsample(wf_list{n}.data, 5);
    f0 = wf_data(2:end,1);
    hpf0 = abs((wf_data(2:end,2) + sqrt(-1).*wf_data(2:end,3)));

    wf_fmin = 1e-5;
    [hpf_fit, f_fit] = FIT_INSPIRAL(wf_fmin, hpf0, f0  );

    % (4a) Epoch-4: LISA

    Detector = 'LISA';
    timeline =4; 


    for i = 1:numel(Mdet)
    
        [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
        [fmin, fmax] = LISA_CBC_FMIN( qq(n), Mdet(i), timeline);
        [rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
        [dAvg_LISA(i,n), zAvg_LISA(i,n)] = AVG_DIST(rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list);
        Msrc_LISA(i,n) = Mdet(i)./(1+zAvg_LISA(i,n));    

    end



    % (4b) Epoch-3: Einstein Telescope

    Detector = 'ET';


    for i = 1:numel(Mdet)
    
        [ hpf, f] = SCALE_WAVEFORM( Mdet(i), Ddet, M0, D0, hpf_fit, f_fit);
        fmin = 1; fmax = 5e3;
        [rho_opt] = SNR_OPTIMAL( hpf, f, Detector, fmin, fmax);
        [dAvg_Ground(i,n), zAvg_Ground(i,n)] = AVG_DIST(rho_opt, Detector, rho_threshold, Ddet, z_list, DL_list);
        Msrc_Ground(i,n) = Mdet(i)./(1+zAvg_Ground(i,n));    

    end
    
    
    
    % (4c) Multiband Distance

    [dMulti(:,n)] = MULTIBAND_DISTANCE(MMulti, dAvg_LISA(:,n), Msrc_LISA(:,n), dAvg_Ground(:,n), Msrc_Ground(:,n));

end


%% (5) INTERPOLATE MULTIBAND-DISTANCE FOR MESH OF MASS-RATIOS & MASSES
clear x y z xyz

z = dMulti(:)./1e3;

x = repmat(MMulti, 1, numel(wf_list));
x = x(:);

y = repmat(eta', 1, numel(MMulti));
y = y(:);
y = sortrows(y, 'descend');

xyz = [x, y, z];
xyz = sortrows(xyz,2);

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

step = 500;
xi=linspace(min(x),max(x),step);
yi=linspace(min(y),max(y),step);

[XI, YI]=meshgrid(xi, yi);

ZI = griddata(x,y,z,XI,YI,'natural');



%% (6) CONTOUR PLOT

figure(1)
numColors = 12;
colormap(flipud(parula())  )
contourf(XI,YI,log10(ZI),'LineStyle','none');

set (gca, 'Xscale', 'log');
axis square
set (gca, 'Yscale', 'log');
set(gca, 'ytick', [0.02; 0.03; 0.05; 0.1; 0.15; 0.25])
caxis([-2, 1.5])
h = colorbar('YTick',linspace(-2,1,4));
xlim([50, 2e4])

xlabel('Total Mass in Source Frame [$M_\odot$]','Interpreter','latex','FontName','TimesNewRoman','FontSize',14), 
ylabel('Symmetric Mass-Ratio','Interpreter','latex','FontName','TimesNewRoman','FontSize',14)




