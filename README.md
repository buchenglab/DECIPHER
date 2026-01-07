DECIPHER is a modular software package for vibrational hyperspectral image enhancement and analysis.

The entire data-processing workflow has been restructured into a streamlined, modular pipeline comprising four customizable stages: hyperspectral image enhancement, background correction, spectral signature identification, and spectral unmixing.

1. Hyperspectral image enhancement. In the first stage, we incorporated methods to correct shading and temporal drift and introduced two denoising algorithms, BM4D and SPEND, for hyperspectral image denoising. These methods mitigate imaging artifacts and improve the signal-to-noise ratio (SNR), enabling more robust spectral unmixing in subsequent steps.
   
2. Background correction across modalities. In the second stage, we implemented modality-specific background removal algorithms for SRS, CARS, and MIP imaging. This step corrects spectral line-shape distortions and improves the signal-to-background ratio (SBR), which is critical for reliable identification of spectral references, particularly in CARS microscopy.

3. Flexible identification of basis spectra. In the third stage, we expanded the framework to include four methods for identifying basis spectra for spectral unmixing. In addition to the original approach guided by pure chemical standards, we introduced three alternative strategies based on sample morphology, spectral phasor analysis, and peak fitting. These additions substantially broaden the range of application scenarios supported by the package.

4. Enhanced spectral unmixing algorithms. In the final stage, the updated DECIPHER package now provides two options for spectral unmixing. Beyond the original pixel-wise LASSO method, and in response to reviewer concerns regarding reliance on predefined basis spectra, we introduced an additional algorithm, MCR-LASSO. This approach enables simultaneous updating of concentration maps and spectral profiles via alternating least-squares optimization, with optional data augmentation to stabilize spectral profile estimation.		

Extensive demonstration datasets are included in the Figshare repository: https://doi.org/10.6084/m9.figshare.31014862
