function [DL,DM,K]=lum_dist(Z,CosmoPars,Lambda,Spectra,SpecType);
%---------------------------------------------------------------------------
% lum_dist function          Compute luminosity distance to an object,
%                          given its redshift and cosmological parameters.
%                          Given the object spectra, calculate also
%                          the K-correction.
% Input  : - Vector of redshifts.
%          - Cosmological parameters : [H0, \Omega_{m}, \Omega_{\Lambda}],
%            default is [70, 0.3, 0.7].
%          - Observation wavelength [Ang]. If SpecType="Hz", then
%            the units are freq. [Hz].
%            If the argument is a string, then it is taken as
%            filter name (check calc_synflux.m for list of filters).
%            Filter name is working only with the 'Ang' option.
%            If this parameter is not given, the K-correction is
%            not calculated.
%          - Object Restframe Spectra, [Wavelength, Specific flux] matrix.
%            If SpecType="Hz", then the units are freq. [Hz].
%            Default is flat power low spectrum: L(nu) = L.
%          - SpecType:
%            'Ang' - Working in wavelength [Ang] - default.
%            'Hz'  - Working in frequency [Hz].
% Output : - Luminosity distance [parsec].
%          - Distance modulus [mag].
%          - K correction [mag].
%            Where, m_{intrinsic} = m_{observed} - K(z)
% Reference : Perlmutter et al. 1997 ApJ, 483, 565
%             Oke & Sandage 1968 ApJ, 154, 21
%             Peterson, B.M., 1997, AGN, p.165
% Tested : Matlab 5.1
%     By : Eran O. Ofek              July 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
C  = 29979245800;    % speed of light [cm/sec]
Pc = 3.0857e18;      % Parsec [cm]

if (nargin==1),
%   CosmoPars = [70, 0.3, 0.7];

   Lambda = NaN;
elseif (nargin==2),
   Lambda = NaN;
elseif (nargin==3),
   % Constant spectra in L(nu)
   Spectra  = [1e-10 1; 1e60 1];
   SpecType = 'Hz'; 
   Lambda   = C./(Lambda*1e-8);
elseif (nargin==4),
   SpecType = 'Ang';
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end


H0     = CosmoPars(1);
OmegaM = CosmoPars(2);
OmegaL = CosmoPars(3);
 

% convert H0 to cm/sec/sec
H0 = H0*100000/(Pc*1e6);



if ((OmegaM+OmegaL)>1),
   Lx = inline('sin(x)','x');
   K  = 1 - OmegaM - OmegaL;
elseif ((OmegaM+OmegaL)<1),
   Lx = inline('sinh(x)','x');
   K  = 1 - OmegaM - OmegaL;
else
   % OmegaM + OmegaL == 1
   Lx = inline('x','x');
   K  = 1;
end

N  = length(Z);
DL = zeros(size(Z));
for I=1:1:N,
   ZI = Z(I);

   Int = inline('((1+ZI).^2.*(1+OmegaM.*ZI) - ZI.*(2+ZI).*OmegaL).^(-1./2)','ZI','OmegaM','OmegaL');

   DL_Int = quad(Int,0,ZI,[],[],OmegaM,OmegaL);

   DL(I) = (C*(1+ZI)/(H0.*sqrt(abs(K)))).*Lx(sqrt(abs(K)).*DL_Int);

end

% luminosity distance in Parsecs:
DL = DL.'./Pc;

% distance modulus:
if (nargout>1),
   DM = 5.*log10(DL./10);
end
DL = DL.';
DM = DM.';

%-------------------------
% K-correction [mag units]
% M = m - K
%-------------------------
if (nargout>2),
   if (isnan(Lambda)==1),
      % do not calculate K-Correction
      K = NaN;
   else
      % calculate K-correction
      Method = 'linear';
      switch SpecType
       case 'Hz'
          % frequency [Hz]
          if (isstr(Lambda)==1),
             error('Filter name is working only in the angs. option');
          else
             Lnu = interp1(Spectra(:,1),Spectra(:,2),Lambda,Method).*ones(size(Z));
             Lnz = interp1(Spectra(:,1),Spectra(:,2),Lambda.*(1+Z),Method);
             K   = -2.5.*log10(1+Z) + 2.5.*log10(Lnu./Lnz);
          end
       case 'Ang'
          % wavelength [Ang]
          if (isstr(Lambda)==0),
             Lla = interp1(Spectra(:,1),Spectra(:,2),Lambda,Method).*ones(size(Z));
             Llz = interp1(Spectra(:,1),Spectra(:,2),Lambda./(1+Z),Method);

             K   = +2.5.*log10(1+Z) + 2.5.*log10(Lla./Llz);
          else
             % integrate using calc_synflux
             K = zeros(N,1);
             for I=1:1:N,
                %--- Note : We should divided the filter transmition
                % by (1+z) instead we multiply the spectrum by (1+z).
                ShiftedSpectra = [Spectra(:,1).*(1+Z(I)), Spectra(:,2)];   %<- doing it in the inverse way
                %ShiftedSpectra = [Spectra(:,1)./(1+Z(I)), Spectra(:,2)./((1+Z(I)).^3)];
                I1   = calc_synflux(Spectra,        0, Lambda, Method);
                I2   = calc_synflux(ShiftedSpectra, 0, Lambda, Method);
                K(I) = +2.5.*log10(1+Z(I)) + 2.5.*log10(I1./I2);
                %K(I) = + 2.5.*log10(I1./I2);
                %K(I) = +2.5.*3.*log10(1+Z(I)) + 2.5.*(I1./I2);
             end
          end
       otherwise
          error('Unknown SpecType option');
      end
      
   end
end



