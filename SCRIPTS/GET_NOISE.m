function [Sn, f_Sn, fmin_Sn, fmax_Sn] = GET_NOISE(Detector)


if (Detector == "LISA")
	curves = importdata('../PSD/L3_noise.txt');
	timeline = 4; %years
	f_Sn = curves(:,1);
	galactic_noise = LISA_GAL_NOISE( curves(:,1), timeline);
	Sn = sqrt(curves(:,2) + galactic_noise);
	fmin_Sn = 2e-5; %Hz
	fmax_Sn = max(f_Sn);

elseif (Detector == "aLIGO")
    curves = importdata('../PSD/Ground_noise.mat');
	f_Sn= curves.aLIGO_design(:,1);
	Sn = curves.aLIGO_design(:,2);
	fmin_Sn = 10;
    fmax_Sn = max(f_Sn);

elseif (Detector == "VIRGO")
    curves = importdata('../PSD/Ground_noise.mat');
    f_Sn= curves.AdVirgo(:,1);
    Sn = curves.AdVirgo(:,2);
    fmin_Sn = 10;
    fmax_Sn = max(f_Sn);


elseif Detector == "KAGRA"
    curves = importdata('../PSD/Ground_noise.mat');
    f_Sn= curves.Kagra(:,1);
    Sn = curves.Kagra(:,2);
    fmin_Sn = 10;
    fmax_Sn = max(f_Sn);

elseif Detector == "A+"
    curves = importdata('../PSD/Ground_noise.mat');
    f_Sn= curves.Aplus(:,1);
    Sn = curves.Aplus(:,2);
    fmin_Sn = 10;
    fmax_Sn = max(f_Sn);

elseif Detector == "VOY"
    curves = importdata('../PSD/Ground_noise.mat');
    f_Sn= curves.Voyager(:,1);
    Sn = curves.Voyager(:,2);
    fmin_Sn = 5;
    fmax_Sn = max(f_Sn);

elseif Detector == "CE"

    curves = importdata('../PSD/Ground_noise.mat');
    f_Sn= curves.CE(:,1);
    Sn = curves.CE(:,2);
    fmin_Sn = 5;
    fmax_Sn = max(f_Sn);


elseif Detector == "ET"
    curves = importdata('../PSD/Ground_noise.mat');
    f_Sn= curves.ET_D(:,1);
    Sn = curves.ET_D(:,2);
    fmin_Sn = 1;
    fmax_Sn = max(f_Sn);

end
end

