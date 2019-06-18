clear , clc


%% (1) GENERATE REDSHIFT VS. LUMINOSITY DISTANCE

CosmoPars = [67.66; 0.3111; 0.6889]; 
% PLANCK 2018 Table 2 (TT, TE, EE +LowE + Lensing + BAO) [H OmegM OmegLamba]

z_list =logspace(-4, log10(5e3), 5000);
[DL_list, ~, ~]  = lum_dist(z_list, CosmoPars );
DL_list = DL_list./1e6; %Get in Mpc


%% (2) LOAD WAVEFORM OF VARRYING MASS-RATIO (qq)

wf_c0   =   importdata('../WF/q-series/f5_q1_M50_D100_i0_IMRHM.dat');                   
wf_cp2  =   importdata('../WF/chi-series/f5_chi0.2_q1_M50_D100_i0_IMRHM.dat');          
wf_cp4  =   importdata('../WF/chi-series/f5_chi0.4_q1_M50_D100_i0_IMRHM.dat');         
wf_cp6  =   importdata('../WF/chi-series/f5_chi0.6_q1_M50_D100_i0_IMRHM.dat');         
wf_cp8  =   importdata('../WF/chi-series/f5_chi0.8_q1_M50_D100_i0_IMRHM.dat');        
wf_cp99 =   importdata('../WF/chi-series/f5_chi0.99_q1_M50_D100_i0_IMRHM.dat');       
wf_cm2  =   importdata('../WF/chi-series/f5_chi-0.2_q1_M50_D100_i0_IMRHM.dat');       
wf_cm4  =   importdata('../WF/chi-series/f5_chi-0.4_q1_M50_D100_i0_IMRHM.dat');       
wf_cm6  =   importdata('../WF/chi-series/f5_chi-0.6_q1_M50_D100_i0_IMRHM.dat');        
wf_cm8  =   importdata('../WF/chi-series/f5_chi-0.8_q1_M50_D100_i0_IMRHM.dat');         
wf_cm99 =   importdata('../WF/chi-series/f5_chi-0.99_q1_M50_D100_i0_IMRHM.dat');        

wf_list = {wf_cp99, wf_cp8, wf_cp6, wf_cp4, wf_cp2, wf_c0, wf_cm2, wf_cm4, wf_cm6, wf_cm8, wf_cm99}
chi =     [0.99,    0.8,    0.6,    0.4,    0.2,    0.0,    -0.2,   -0.4,   -0.6,   -0.8,   -0.99];


qq=1;
eta = qq./(1+qq).^2;
M0 = 50; % wf total-mass in Msun
D0 = 100; % wf distance in Mpc

%% (3) CREATE VARIABLES

Mdet=logspace(2,8, 1000); % from 10^1 to 10^7 solar mass
Ddet=100; % Mpc
rho_threshold  = 8;

Msrc_LISA = zeros(numel(Mdet),numel(wf_list));
dAvg_LISA = zeros(numel(Mdet),numel(wf_list));
zAvg_LISA = zeros(numel(Mdet),numel(wf_list));

Msrc_Ground = zeros(numel(Mdet),numel(wf_list));
dAvg_Ground = zeros(numel(Mdet),numel(wf_list));
zAvg_Ground = zeros(numel(Mdet),numel(wf_list));

MMulti = logspace(log10(1e3), log10(2e4), 3000); % from 50 to 2*10^5 solar mass
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
        [fmin, fmax] = LISA_CBC_FMIN( qq, Mdet(i), timeline);
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
    
end

for n=1:numel(wf_list)    
    % (4c) Multiband Distance

    [dMulti(:,n)] = MULTIBAND_DISTANCE(MMulti, dAvg_LISA(:,n), Msrc_LISA(:,n), dAvg_Ground(:,n), Msrc_Ground(:,n));

end


%% (5) INTERPOLATE MULTIBAND-DISTANCE FOR MESH OF MASS-RATIOS & MASSES
clear x y z xyz ZI XI YI

z = dMulti(:)./1e3;

x = repmat(MMulti, 1, numel(wf_list));
x = x(:);

y = repmat(chi', 1, numel(MMulti));
y = y(:);
y = sortrows(y, 'descend');

xyz = [x, y, z];
xyz = sortrows(xyz,2);

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

step = 3000;
xi=logspace(log10(min(x)),log10(max(x)),step);
yi=linspace(min(y),max(y),step);

[XI, YI]=meshgrid(xi, yi);

ZI = griddata(x,y,z,XI,YI,'natural');



%% (6) CONTOUR PLOT

figure(1)
numCont = 9;
colormap(flipud(parula())  )
contourf(XI,YI,log10(ZI),numCont,'LineStyle','none');

set (gca, 'Xscale', 'log');
axis square
%set(gca, 'ytick', [-0.8; -0.6; -0.4; -0.2; 0.0; 0.2; 0.4; 0.6; 0.8])

cmin = -2;
cmax = 1.5;
caxis([cmin, cmax])
h = colorbar('YTick',linspace(cmin,cmax, 8));
xlim([1e3, 2e4])

xlabel('Total Mass in Source Frame [$M_\odot$]','Interpreter','latex','FontName','TimesNewRoman','FontSize',14), 
ylabel('Effective Spin','Interpreter','latex','FontName','TimesNewRoman','FontSize',14)




