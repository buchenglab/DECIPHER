clc;
clear;

tic;
%%
Spec_filename = 'Single_Spec_sample.txt';
Spec = load(Spec_filename);

Spec = Spec(:,2)';

bkg_ar = arPLS_baseline_v0(Spec,1e2,1e-6)';
peak_ar = Spec - bkg_ar;
peak_ar_norm = (peak_ar-min(peak_ar))./(max(peak_ar) - min(peak_ar));

figure;
plot(bkg_ar,'DisplayName','bkg');hold on
plot(peak_ar,'DisplayName','signal');
plot(Spec,'DisplayName','Raw Spec');
hold off;
legend();

figure;
plot(peak_ar_norm);