clc;
clear;

tic;
%%
file_raw = 'HEK293T H2B GFP_100ns_g1kmW_600k_20x__FOV 0_250X250X200nm_dwell40us_16h25m13s_Hyper_AutoWN0_LIA.tif';
img_raw = double(tiffreadVolume(file_raw));
[xx, yy, l] = size(img_raw);


%%
[img_raw_arPLS_txt, bkg_raw_arPLS_txt] = arPLS_process(img_raw,1e2,1e-6,1);
dlmwrite([file_raw 'raw_arPLS_test.txt'],img_raw_arPLS_txt);
dlmwrite([file_raw 'raw_arPLS_test_bkg.txt'],bkg_raw_arPLS_txt);



toc;

%%
function [img_arPLS, bkg_arPLS] = arPLS_process(img,val_smooth,val_diff,n)
    [xx,yy,l] = size(img);
    img_arPLS = zeros([xx yy l]);
    bkg_arPLS = zeros([xx yy l]);
    parfor ii=1:xx
        for jj=1:yy
            img_pro = img(max(ii-n,1):min(xx,ii+n),max(jj-n,1):min(yy,jj+n),:);
            spec = mean(img_pro,[1 2]);
            spec = spec(:,:);
            raw = img(ii,jj,:);
            raw = raw(:,:);
            bkg_ar = arPLS_baseline_v0(spec,val_smooth,val_diff)';
            peak_ar = raw - bkg_ar;
            img_arPLS(ii,jj,:) = peak_ar;
            bkg_arPLS(ii,jj,:) = bkg_ar;
        end
    end
end
%%
function [img_airPLS, bkg_airPLS] = airPLS_process(img,val_smooth,val_diff,n)
    [xx,yy,l] = size(img);
    img_airPLS = zeros([xx yy l]);
    bkg_airPLS = zeros([xx yy l]);
    for ii=1:xx
        for jj=1:yy
            img_pro = img(max(ii-n,1):min(400,ii+n),max(jj-n,1):min(400,jj+n),:);
            spec = mean(img_pro,[1 2]);
            spec = spec(:,:);
            raw = img(ii,jj,:);
            raw = raw(:,:);
            [peak_air,bkg_air] = airPLS(spec,1e3,3,0.9,0.1,50);
            peak_air = raw-bkg_air;
            img_airPLS(ii,jj,:) = peak_air;
            bkg_airPLS(ii,jj,:) = bkg_air;
        end
    end
end