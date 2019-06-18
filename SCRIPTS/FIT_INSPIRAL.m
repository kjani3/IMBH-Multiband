function [ hpf_fit, f_fit] = FIT_INSPIRAL(fmin, hpf, f  )
% FIT_INSPIRAL fits the freqeuncy-domain inspiral signal to a given fmin

stp=100;
[~,I] = max(hpf);
x = f(I:I+stp);
y = hpf(I:I+stp);
xx = logspace(log10(fmin), log10(x(end)), 1e6);
[slope, intercept] = logfit(x,y,'loglog');
yApprox = (10^intercept)*xx.^(slope);

f_fit = [xx'; f(I+stp+1:end)];
hpf_fit =[yApprox';  hpf(I+stp+1:end)];

end

