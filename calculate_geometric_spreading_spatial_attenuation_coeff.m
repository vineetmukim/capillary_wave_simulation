clc; clear all; close all;
tic;
freq = 20; % oscillation freq
omega = 2*pi*freq; % angular freq
time_list = 0:1e-3:1; % simulated time
x_list = (0.1:1e-2:80)*1e-3; % radius
C = 1e-4; % amp
lambda = 10e-3; % wavelength
k = 2*pi/lambda;
attenuation_coeff_list = nan(length(time_list), 1);

% fig = figure();
% ax = axes(fig);
parfor i=1:length(time_list)
    disp(time_list(i));
%     cla(ax);
    y_list = C*sin(k*x_list-omega*time_list(i))*1./sqrt(x_list);

%     plot(ax, x_list, y_list);
%     hold on;
%     ylim([-20*C, 20*C]);

    [~, peak_indices] = findpeaks(abs(y_list));
    x_peaks = x_list(peak_indices);
    y_peaks = y_list(peak_indices);
%     scatter(ax, x_peaks, y_peaks, 'r*');
%     x_peaks_norm = x_peaks-x_peaks(1);
    % fit
    [fitresult, ~] = fit(x_peaks', abs(y_peaks'), 'exp1');
    attenuation_coeff_list(i) = fitresult.b;
    
%     pause(0.05);
end
toc;
figure;
scatter(time_list, attenuation_coeff_list);
disp(['mean attenuation coeff = ', num2str(mean(attenuation_coeff_list))]);
