function [cbc_fmin, cbc_fmax] = LISA_CBC_FMIN( q, M, timeline)
% LISA_CBC_FMIN computes the frequency spanned by a binary source during
% the mission lifetime of LISA
%   

SEC_TO_YR = 31536000;
G = 6.67408e-11;
c = 3e8;
Msun = 2e30;

eta = q./(1+q).^2;
Mchirp = eta.^(3/5).*M.*Msun;
f1 = 1/(8*pi*G*Mchirp/c^3);
f2 = (5*G*Mchirp/c^3)^(3/8);
f3min = 1/(timeline.*SEC_TO_YR)^(3/8);
cbc_fmin = f1*f2*f3min;
cbc_fmax = 1;


if cbc_fmin<2e-5
    cbc_fmin = 2e-5;
end



end



