function [ dAvg, zAvg] = AVG_DIST(rho_opt, Detector, rho_th, d0, z_list, DL_list)
% AVG_DIST computes sky-angle and inclination average distance at S/N threshold rho_th in a Detector.
% rho_opt refers to the S/N measured for a source at distance d0 Mpc

if (Detector == "LISA")

	angle_avg_factor = sqrt(4/5);

else
	angle_avg_factor = sqrt(4/25);
end


dAvg =angle_avg_factor.*d0.*(rho_opt./rho_th);       

[~,I] = min(abs(dAvg-DL_list));

zAvg = z_list(I);    

end

