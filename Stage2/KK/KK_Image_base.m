clc;
close all
clear all;
%%
tic;
file_name = '2_my_model_20250907_OPA_CARS_FP_T24_4unet_100uM.tif';
img = double(tiffreadVolume(file_name));
img = img(:,:,1:99);

[xx,yy,m]=size(img);
Total_L = m;
img_kk = zeros(xx,yy,m);
img_kk_arPLS = zeros(xx,yy,m);
img_kk_asls = zeros(xx,yy,m);

REF_filename = "BKG_rich_2avg.txt";
I_REF = load(REF_filename);
I_REF = I_REF(1:m,2);
I_REF = (I_REF-min(I_REF))/(max(I_REF)- min(I_REF))+2;


parfor ii = 1:xx
    for jj = 1:yy
        I_CARS = squeeze(img(ii,jj,:));

        I_CARS = (I_CARS-min(I_CARS))/(max(I_CARS)- min(I_CARS))+2;

        [KK_spec_imag, KK_spec_real] = KKHilbert(I_REF, I_CARS);
        KK_ref = KK_spec_real + 1j*KK_spec_imag;

        phase_ref = angle(KK_ref);
        
        smooth_fac = 5e1;

        phase_error_arpls = arPLS_baseline_v0(phase_ref, 1e4, 1e-6);
        amp_corr_arpls = exp(imag(Hilbert(phase_error_arpls)));
        scale_arpls = 1./sgolayfilt(amp_corr_arpls.*abs(KK_ref).*cos(phase_ref-phase_error_arpls),3,Total_L);
        spec_kk_arPLS = amp_corr_arpls.*scale_arpls.*abs(KK_ref).*sin(phase_ref - phase_error_arpls);
        img_kk_arPLS(ii,jj,:) = spec_kk_arPLS;

        phase_error_asls = asLS_baseline_v1(phase_ref, 1e4, 1e-4);
        amp_corr_asls = exp(imag(Hilbert(phase_error_asls)));
        scale_asls = 1./sgolayfilt(amp_corr_asls.*abs(KK_ref).*cos(phase_ref-phase_error_asls),3,Total_L);
        spec_kk_asls = amp_corr_asls.*scale_asls.*abs(KK_ref).*sin(phase_ref - phase_error_asls);
        img_kk_asls(ii,jj,:) = spec_kk_asls;
        ii
    end
end

%%
dlmwrite(['2_Sub_KK_PLS_' file_name '.txt'], img_kk_arPLS);
dlmwrite(['2_Sub_KK_AsLS_' file_name '.txt'], img_kk_asls);
toc;
%%
function peak_arPLS = arPLS_spec(raw, smooth_v, diff_v)
    bkg_arPLS = arPLS_baseline_v0(raw,smooth_v,diff_v)';
    peak_arPLS = raw - bkg_arPLS;
end