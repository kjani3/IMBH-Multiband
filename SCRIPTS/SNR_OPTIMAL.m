function [ rho_opt] = SNR_OPTIMAL( wf, f_wf, Detector, fmin, fmax)
% SNR_OPTIMAL computes the optimal S/N in a Detector for a waveform wf with
% GW frequnecy f_wf. The integral is from fmin to fmax


 [Sn, f_Sn, ~, ~] = GET_NOISE(Detector);
 
 [~, I_fmin]   =  min(abs(f_wf - fmin));
 I_fmin = I_fmin +1;
 [~, I_fmax]    =  min(abs(f_wf - fmax));
 I_fmax = I_fmax - 1;
 
 flim =  f_wf(I_fmin:I_fmax);
 wflim = wf(I_fmin:I_fmax);

 Snlim = interp1(f_Sn, Sn.^2, flim,'spline');
 func_int = wflim.*conj(wflim)./Snlim;
 

 rho_opt= sqrt( 4.*trapz(flim, func_int ) );
 
 
 
end



