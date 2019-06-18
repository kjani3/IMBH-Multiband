function [ hpf_scale, f_scale] = SCALE_WAVEFORM( M_scale, D_scale, M0, D0, hpf0, f0)
% SCALEE WAVEFORM scales a frequency domain IMR waveform generated at 
% [mass, distance] =[M0, D0] to  [M_scale, D_scale]


freq_scale = M0./M_scale;
amp_scale = (M_scale./M0).^2.*(D0./D_scale);

f_scale = f0.*freq_scale;
hpf_scale = hpf0.*amp_scale;

end

