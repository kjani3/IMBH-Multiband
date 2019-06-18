function [ dHor, zHor] = HOR_DIST( rho_opt, Detector, rho_th, d0, z_list, DL_list )
% HOR_DIST computes horizon distance at S/N threshold rho_th in a Detector.
% rho_opt refers to the S/N measured for a source at distance d0 Mpc


if (Detector == "LISA")
    opt_factor = sqrt(5);
else
    opt_factor = 1;
end
    
    
dHor =opt_factor.*d0.*(rho_opt./rho_th);    

[~,I] = min(abs(dHor-DL_list));

zHor = z_list(I);    

end

