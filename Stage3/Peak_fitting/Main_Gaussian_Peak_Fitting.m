%% Gaussian Peak Fitting (MATLAB version)
% Requires: Curve Fitting Toolbox

close all; clear; clc;

%% Load and preprocess data (match Python logic)

y = load('Raw_Spec_1600.txt');   % assumes numeric text file with >=2 cols

% Normalize: (y - min) / (max - min)
y = (y - min(y)) ./ (max(y) - min(y));
x = linspace(1,size(y,1),size(y,1))';

%% Define Gaussian model
% 3-Gaussian model is provided as an example. Users should change accordingly.
% lmfit GaussianModel uses:
% y = amplitude/(sigma*sqrt(2*pi)) * exp(-(x-center)^2/(2*sigma^2))
%

gauss = @(A, c, s, x) (A./(s*sqrt(2*pi))) .* exp(-((x - c).^2)./(2*s.^2));

% Specify the number of peaks, and peak amplitude, center, and width
ft = fittype(@(g1_amplitude,g1_center,g1_sigma, ...
               g2_amplitude,g2_center,g2_sigma, ...
               g3_amplitude,g3_center,g3_sigma, x) ...
    gauss(g1_amplitude,g1_center,g1_sigma,x) + ...
    gauss(g2_amplitude,g2_center,g2_sigma,x) + ...
    gauss(g3_amplitude,g3_center,g3_sigma,x), ...
    'independent', 'x', ...
    'coefficients', {'g1_amplitude','g1_center','g1_sigma', ...
                     'g2_amplitude','g2_center','g2_sigma', ...
                     'g3_amplitude','g3_center','g3_sigma'});

% Initial guess of peak position and width (from your Python code)
startVals = [ ...
    0.1, 26, 1, ...   % g1
    0.5, 51, 1, ...   % g2
    0.6, 40, 1  ...   % g3 
];

opts = fitoptions(ft);
opts.StartPoint = startVals;

% Important: add bounds for stability. Peak position bound is recommended
 opts.Lower = [0, 20, 0,  ...
               0, 40,  0,  ...
               0, 30, 0  ];

 opts.Upper = [Inf, 30, Inf, ...
               Inf, 60, Inf, ...
               Inf, 50, Inf];

%% Fit
[fitObj, gof] = fit(x, y, ft, opts);

disp('=== Fit coefficients ===');
disp(fitObj);
disp('=== Goodness of fit ===');
disp(gof);

% Create fitted curve
y_fit = feval(fitObj, x);

%% Evaluate individual components (like result.eval_components())
c = coeffvalues(fitObj);

g1 = gauss(c(1), c(2), c(3), x);
g2 = gauss(c(4), c(5), c(6), x);
g3 = gauss(c(7), c(8), c(9), x);

%% Save components (matches your Python naming style)
writematrix(g1, 'g1_.txt', 'Delimiter', 'tab');
writematrix(g2, 'g2_.txt', 'Delimiter', 'tab');
writematrix(g3, 'g3_.txt', 'Delimiter', 'tab');

%% Plot
figure('Position', [200 200 900 500]);
plot(x, y, 'b-', 'LineWidth', 1.5); hold on;
plot(x, y_fit, 'r-', 'LineWidth', 1.5);

plot(x, g1, '--', 'LineWidth', 1.2);
plot(x, g2, '--', 'LineWidth', 1.2);
plot(x, g3, '--', 'LineWidth', 1.2);

legend({'Observed Spectrum','Fitted Spectrum','g1_ (Gaussian)','g2_ (Gaussian)','g3_ (Gaussian)'}, ...
       'Location','best');
xlabel('x');
ylabel('Intensity');
title('Gaussian Peak Fitting (MATLAB)');
grid on;
