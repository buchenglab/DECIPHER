clc;
close all
clear;
%%


WN = linspace(1476.75,1930,100)-20;
WN = linspace(21,120,200)-20;
% I_REF = load("BKG_rich.txt");

REF_filename = "BKG_rich.txt";

I_REF = load(REF_filename);
I_REF = I_REF(1:199,2);
I_REF = (I_REF-min(I_REF))/(max(I_REF)- min(I_REF))+3;
% I_REF = (I_REF-I_REF(1))/(I_REF(150)- I_REF(1))+3;

% filename = 'BSA.txt';
% I_CARS= load(filename);

CARS_filename = 'BSA.txt';
I_CARS= load(CARS_filename);

I_CARS = I_CARS(1:199,2);
I_CARS = (I_CARS-min(I_CARS))/(max(I_CARS)- min(I_CARS))+3;
% I_CARS = (I_CARS-I_CARS(1))/(I_CARS(150)- I_CARS(1))+3;

Total_L = length(I_CARS);
WN = WN(1:Total_L);
I_NRB = I_REF;
%% Basic Retrieval

[KK_ideal_imag, KK_ideal_real] = KKHilbert(I_NRB, I_CARS);
KK_ideal = KK_ideal_real + 1j*KK_ideal_imag;
phase_ideal = angle(KK_ideal);

figure
plot(WN,I_CARS, 'LineWidth', 2);
hold all
plot(WN, I_NRB,'LineWidth',2);
legend('CARS','NRB');
xlabel('Wavenumber (cm^{-1})');
ylabel('Signal Int (au)');
title('CARS and NRB Signal (Spectra)');

figure
plot(WN, KK_ideal_imag,'LineWidth',2);
hold all
legend('Retrieved');
xlabel('Wavenumber (cm^{-1})');
ylabel('Raman-Like Int. (no units)');
title('Raman-Like Spectra (Retrieved)');

% outname = ['Spec_' filename];
% dlmwrite(outname,KK_ideal_imag);

%%
[KK_ref_imag, KK_ref_real] = KKHilbert(I_REF, I_CARS);

KK_ref = KK_ref_real + 1j*KK_ref_imag;
phase_ref = angle(KK_ref);

figure
plot(WN,I_REF)
hold all
plot(WN, I_NRB+0.1)
legend('Surrogate NRB', 'NRB');
xlabel('Wavenumber (cm^{-1})');
ylabel('Signal Int. (au)');
title('Ideal and Surrogate NRB Spectra');


figure
% plot(WN,KK_ref_imag,'LineWidth',2);
hold all
plot(WN,phase_ref)
legend('Retrieved ')
xlabel('Wavenumber (cm^{-1})');
ylabel('Raman-Like Int. (no units)');
title('Raman-Like Spectra with Surrogate NRB Spectrum');

%% Phase-error correction (Note: 'ideal' now means ideal correction, 
%  not the actual ideal [i.e., without any error]
phase_error_asls = asLS_baseline_v1(phase_ref, 1e2, 1e-4);
phase_error_arpls = arPLS_baseline_v0(phase_ref, 1e3, 1e-6);

figure
subplot(2,1,1)
plot(WN, phase_ref);
hold all
plot(WN, phase_ideal);
legend('Phase (Ref)', 'Phase (NRB, Ideal)')
xlabel('Wavenumber (cm^{-1})');
ylabel('Phase (rad)');
title('Phase');

subplot(2,1,2)
plot(WN,phase_error_asls)
hold all
plot(WN,phase_error_arpls)
legend('Phase error (AsLS)', 'Phase error (arPLS)')
xlabel('Wavenumber (cm^{-1})');
ylabel('Phase (rad)');
title('Phase Error');

% amplitude correction factors
amp_corr_asls = exp(imag(Hilbert(phase_error_asls)));
amp_corr_arpls = exp(imag(Hilbert(phase_error_arpls)));

%% Scaling and final plots
% scale_ideal = (KK_ideal_real./(amp_corr_ideal.*abs(KK_ref).*cos(phase_ref-phase_error_ideal)));
scale_asls = 1./sgolayfilt(amp_corr_asls.*abs(KK_ref).*cos(phase_ref-phase_error_asls),3,Total_L);
scale_arpls = 1./sgolayfilt(amp_corr_arpls.*abs(KK_ref).*cos(phase_ref-phase_error_arpls),3,Total_L);
% 
figure
plot(WN, amp_corr_asls.*scale_asls.*abs(KK_ref).*sin(phase_ref - phase_error_asls),'LineWidth',2)
hold all
plot(WN, amp_corr_arpls.*scale_arpls.*abs(KK_ref).*sin(phase_ref - phase_error_arpls))
% plot(WN, amp_corr_ideal.*scale_ideal.*abs(KK_ref).*sin(phase_ref - phase_error_ideal),'r')
% % plot(WN, imag(CHI_R)./abs(CHI_NR),'k');
% 
legend('Correction', 'Correction (arPLS)','Ideal Correction', 'Ideal')
xlabel('Wavenumber (cm^{-1})');
ylabel('Raman-Like Int. (no units)');
title('Raman-Like Spectra with Phase and Amplitude Error Correction');
% xlim([2050,2220]);

Spec_asls = amp_corr_asls.*scale_asls.*abs(KK_ref).*sin(phase_ref - phase_error_asls);
outname_asls = ['AsLS_Spec_' CARS_filename];
dlmwrite(outname_asls,Spec_asls);

Spec_arpls = amp_corr_arpls.*scale_arpls.*abs(KK_ref).*sin(phase_ref - phase_error_arpls);
outname_arpls = ['arPLS_Spec_' CARS_filename];
dlmwrite(outname_arpls,Spec_arpls);
