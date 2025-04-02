clear; clc;
close all;
tic;
% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');

% name_list = ["pl\f\10", "pl\f\20", "pl\f\30", "pl\f\40", "pl\f\50",...
%     "pl\sigma\0.02_20", "pl\sigma\0.10_20","pl\sigma\0.20_20", "pl\sigma\0.50_20",...
%     "pl\mu\10e-3_20", "pl\mu\25e-3_20", "pl\mu\50e-3_20", "pl\mu\100e-3_20",...
%     "pl\mu\10e-3_40", "pl\mu\25e-3_40", "pl\mu\50e-3_40", "pl\mu\100e-3_40"];

% name_list = ["ax\f\10", "ax\f\20", "ax\f\30", "ax\f\40", "ax\f\50",...
%     "ax\sigma\0.02_20", "ax\sigma\0.10_20","ax\sigma\0.20_20", "ax\sigma\0.50_20",...
%     "ax\mu\10e-3_20", "ax\mu\25e-3_20", "ax\mu\50e-3_20", "ax\mu\100e-3_20",...
%     "ax\mu\10e-3_40", "ax\mu\25e-3_40", "ax\mu\50e-3_40", "ax\mu\100e-3_40"];

% name_list = ["pl\mu\10e-3_20"];%, "pl\mu\25e-3_20", "pl\mu\50e-3_20", "pl\mu\100e-3_20"];%, "pl\mu\100e-3_40"];%, "pl\mu\10e-3_20", "pl\mu\10e-3_25", "pl\mu\50e-3_20", "pl\mu\100e-3_20"];


top_num = 9; % top n+1 peaks i.e. temporal indices for averaging
deltaT = 0.001;
folder_path = ['C:\Users\mukim\Documents\PhD\Work\Projects\Capillary Wave Surface Reconstruction\Journal Paper 3\pl_vs_ax_dataset\'];
tic

syms omega k nu sigma;

rho = 1e3; % kg/m^3
g = 9.81; % m/s^2
freq = 10; % Hz
sigma_val = 0.07; % N/m
mu_val = 1e-3; % Pa.s

% full dispersion equation
F(omega, k, nu, sigma) = ((-1i*omega+2*nu*k^2)^2+g*k+sigma*k^3/rho)^2-(4*nu^2*k^3)^2*(k^2-1i*omega/nu); % sq terms


for par_id=1:length(name_list)
    parameter = name_list(par_id);
    data_path = strcat(folder_path, parameter);
    saving_path = strcat(data_path, '\processed_data\');
    
    if contains(parameter, '\f\')
        id_temp = strfind(parameter, '\');
        p = char(parameter);
        freq = str2double(p(id_temp(end)+1:end));
    end
    
    if contains(parameter, '\sigma\')
        id_temp1 = strfind(parameter, '\');
        id_temp2 = strfind(parameter, '_');
        p = char(parameter);
        sigma_val = str2double(p(id_temp1(end)+1:id_temp2-1));
        freq = str2double(p(id_temp2+1:end));
    end   
    
    if contains(parameter, '\mu\')
        id_temp1 = strfind(parameter, '\');
        id_temp2 = strfind(parameter, '_');
        p = char(parameter);
        mu_val = str2double(p(id_temp1(end)+1:id_temp2-1));
        freq = str2double(p(id_temp2+1:end));
    end  
    

    if ~isfolder(saving_path)
        disp('processed folder not found...creating the folder');
        mkdir(saving_path);
    end

    % write plot in video file
    % par_index = strfind(parameter, '\');
    % name = parameter(par_index(1)+1:end);
    name = strrep(parameter,'\','_');

    axisymm_flag = contains(parameter, 'ax\');
    mu_flag = contains(parameter, '\mu\');
    
    files = dir(saving_path);
    
    for i=1:length(files)
    %     disp(i);
    %     disp(files(i).name);
        if contains(files(i).name, 'time') || contains(files(i).name, 'mu') || contains(files(i).name, 'space')
            disp('deleting existing post-processed files');
            %         disp(files(i).name);
            delete(strcat(files(i).folder, '\', files(i).name));
        end
    end
    
    time_end = 0;  
    data_files = dir(data_path);
    for i=1:length(data_files)    
%         disp(data_files(i).name);
        if ~isnan(str2double(data_files(i).name)) && str2double(data_files(i).name)<=1  
            time_end = max(str2double(data_files(i).name), time_end);
        end
    end
%     disp(time_end);
    time_end = floor(time_end*10)/10;
    time_steps = time_end-0.1:deltaT:time_end;
    
    peak_num=0; % for all other cases
 
%     if contains(name, '100e-3')
%         peak_num = 7;
%     elseif contains(name, '50e-3')
%         peak_num = 10;
%     elseif contains(name, '25e-3')
%         peak_num = 10;
%     end
    %-------------------------------------------------------------
    % plot figure size and location
    fig = figure;
    fig.Position = [30 30 1200 400];
    ax = axes(fig);  
    set(gcf, 'Color', [1 1 1]); % set white background to fig    
    pbaspect([4,1,1]);

    space_max_amp_list = nan(length(time_steps), 1);
    space_avg_lambda_list = nan(length(time_steps), 1);
    space_avg_k_list = nan(length(time_steps), 1);    
    space_attenuation_coeff_list = nan(length(time_steps), 2);   

    for i=1:length(time_steps)    
        time = time_steps(i);
        disp(strcat(parameter, '---', num2str(time), ' s'));
        data = [];
        time_file = strcat(folder_path, parameter, '\', num2str(time), '\freeSurf\points');
    %     disp(time_file);
        fid = fopen(time_file, 'r');
        tline = fgetl(fid);
        while ischar(tline)
            num = sscanf(tline, '(%f %f %f)');
    %         disp(size(num));
            if(length(num)==3 && num(3)~=1e-05)
    %             disp(num);
                data = [num data];
            end
            tline = fgetl(fid);
        end
        fclose('all');

        x = data(1, :);
        y = data(2, :)-0.01; % making undisturbed level 0 

        cla(ax);
        scatter(ax, x, y, '.');
        
        
        if axisymm_flag
            xlabel(ax, 'radius (m)');
        else
            xlabel(ax, 'length (m)');   
        end
        ylabel(ax, 'amplitude (m)');
        pbaspect([4 1 1]);

        xlim(ax,[0 80*1e-3]);
        if axisymm_flag
            y_lim_range = [-0.05 0.05]*1e-3;  
            ax.YAxis.TickValues = y_lim_range(1):mean(abs(y_lim_range)):y_lim_range(2);
        else
            y_lim_range = [-0.2 0.2]*1e-3; 
            ax.YAxis.TickValues = y_lim_range(1):mean(abs(y_lim_range)):y_lim_range(2);            
        end
        ylim(ax, y_lim_range);
        grid (ax, 'on');
        set(ax,'XMinorGrid', 'on', 'YMinorGrid', 'on');

        %---------------------------------------------------

        % fitting spline
        coordinates = sortrows([x;y]');
        opts = fitoptions('Method', 'SmoothingSpline'); 
        opts.SmoothingParam = 0.999999999; %0.9999999999;  
        [fitresult_space_avg_amp, ~] = fit(coordinates(:, 1), coordinates(:, 2), fittype('smoothingspline'), opts);   
        x_fit = [1.5:1e-2:80.0]*1e-3;
        y_fit = fitresult_space_avg_amp(x_fit');
        hold (ax, 'on');
        plot(ax, x_fit, y_fit, 'r.');
        title(['time = ', num2str(time*1000), ' ms'], 'Fontweight', 'normal');
        set(ax, 'FontSize', 12, 'FontName', 'Times');
        
        %---------------------------------------------------
        % finding peaks, number of wavelengths and spatial attenuation

%         peak_id_positive_y = peak_id_all(find(y_fit(peak_id_all)>=0));
%         [~, peak_id_positive_y] = findpeaks(y_fit, 'MinPeakHeight', 0); % ensuring positive y
        [~, peak_id_positive_y] = findpeaks(y_fit, 'MinPeakHeight', 0.02*max(abs(y_fit))); % ensuring positive y
        
        [~, peak_id_all] = findpeaks(abs(y_fit));%, 'MinPeakHeight', 0.05*max(abs(y_fit)));
        
        if contains(name, 'ax')
            peak_id_all = [1; peak_id_all];
        end
        
        % limiting number of peaks for cases where peaks with y very close
        % to 0 can change attenuation coeff...mu 1e-1 and 5e-2 and axisymm
        % sigma 0.50 cases
        if peak_num~=0 && length(peak_id_all)>peak_num 
            peak_id_all = peak_id_all(1:peak_num);
            if round(peak_num/2)<=length(peak_id_positive_y)
                peak_id_positive_y = peak_id_positive_y(1:round(peak_num/2));
            end
        end


        
        
        
        
        peak_id_list = peak_id_positive_y;
%         peak_id_list = peak_id_all; % no need to sort data since prepareDataCurve used
        % before calculating the attenuation coeff
%         peak_id_all_with_missing = sort([peak_id_all; missing_id]);      
        

        peak_x_list = (x_fit(peak_id_list))'; % column vector
        peak_y_list = (y_fit(peak_id_list));
        
        peak_x_list_lambda = x_fit(peak_id_positive_y)'; % column vector for lambda
        peak_y_list_lambda = y_fit(peak_id_positive_y)';

        if length(peak_x_list_lambda)>=2
            % space avg lambda in mm
            space_avg_lambda_list(i) = (peak_x_list_lambda(end)-peak_x_list_lambda(1))*1e3/(length(peak_x_list_lambda)-1);     
            space_avg_k_list(i) = 2*pi*1e3/space_avg_lambda_list(i);
        end
        %---------------------------------------------------
        hold(ax, 'on');
        scatter(ax, peak_x_list, peak_y_list, 'ko');
        
%         if mu_flag==1
%             scatter(ax, peak_x_list, peak_y_list, 'b*');  
%         else
%             scatter(ax, peak_x_list_lambda, peak_y_list_lambda, 'b*');   
%         end

        %---------------------------------------------------
        % fitting and finding attenuation coeff
        
        if length(peak_id_list)>=2
%         disp('spatial attenuation fit');

%             space_peak_max_amp_list(i) = max(abs(peak_y_list(find(peak_x_list<0.02)))); 
%             space_peak_max_amp_list(i) = max(abs(peak_y_list)); 
            space_max_amp_list(i) = max(abs(y_fit)); % for axisymm mu 1e-1 case          
            
            % sorting peaks data
            [peak_x_list, peak_y_list] = prepareCurveData(peak_x_list, peak_y_list);               
      
            % normalizing x data
%             peak_x_list_norm = peak_x_list - peak_x_list(1);    
            peak_x_list_norm = peak_x_list; % old way i.e. without normalizing
            
            [fitresult_attenuation_spatial, ~] = fit(peak_x_list_norm, abs(peak_y_list), 'exp1'); 
%             disp(fitresult_attenuation_spatial.a); disp(fitresult_attenuation_spatial.b);
            space_attenuation_coeff_list(i, 1) = fitresult_attenuation_spatial.a;
            space_attenuation_coeff_list(i, 2) = fitresult_attenuation_spatial.b;
            
%             scatter(ax, peak_x_list_norm+peak_x_list(1), abs(peak_y_list));
            if space_attenuation_coeff_list(i, 2)>=0
                % have to change sign of attenuation coeff to match with
                % e^(ikx-wt) formulation, multiplication by i changes the
                % sign of attenuation coeff that actually results in
                % amplitude attenuation and not amplification
                t1 = sprintf('k = %0.1f - %0.1f i', space_avg_k_list(i), abs(space_attenuation_coeff_list(i, 2)));
            else
                t1 = sprintf('k = %0.1f + %0.1f i', space_avg_k_list(i), abs(space_attenuation_coeff_list(i, 2)));
            end 
%             if mu_flag==1
%                 disp('mu flag found...showing amp attenuation line');
                % in fit below, modified x axis is used to correct for
                % normalization axis used in finding the attenuation coeff
%                 y_fit_spatial =
%                 fitresult_attenuation_spatial(x_fit-peak_x_list(1)); %
%                 for normalised x axis
                y_fit_spatial = fitresult_attenuation_spatial(x_fit);
                plot(ax, x_fit, y_fit_spatial, 'k--');%, 'LineWidth', 1.5);
                text(max(x_fit)*0.45, y_lim_range(2)*0.9, sprintf('space averaged wavelength = %0.1f mm', space_avg_lambda_list(i)), 'FontName','Times', 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
%                 text(max(x_fit)*0.55, y_lim_range(2)*0.9, sprintf('y=%0.1de^{%0.1fx}', space_attenuation_coeff_list(i, 1), space_attenuation_coeff_list(i, 2)), 'FontName','Times', 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle'); 
                text(max(x_fit)*0.55, y_lim_range(2)*0.9, t1, 'FontName','Times', 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');                
                legend(ax, {'interface (OpenFOAM)', 'interface (smoothened)', 'crest or trough', 'spatial attenuation fit'}, 'Location', 'south', 'Orientation', 'horizontal');
%             else
%                 % for param sweeps with f and sigma
%                 text(max(x_fit)*0.5, y_lim_range(2)*0.9, sprintf('space averaged wavelength = %0.1f mm', space_avg_lambda_list(i)), 'FontName','Times', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');        
%                 legend(ax, {'interface (OpenFOAM)', 'interface (smoothened)', 'crest or trough'}, 'Location', 'south', 'Orientation', 'horizontal');
%             end
            t1 = t1(find(~isspace(t1)));
            save_title = ['time=', num2str(time*1000), 'ms, max amp=', num2str(space_max_amp_list(i)*1e3,'%0.4f'), 'mm, space avg wavelength=', num2str(round(space_avg_lambda_list(i), 2), '%0.1f'), 'mm, ', t1];
        else
            % exceptional cases with less than 2 points in peak_id_all
            save_title = ['time=', num2str(time*1000), 'ms'];%, space peak max amp=', num2str(round(space_max_amp_list(i)*1e3,4)), 'mm, space avg wavelength=', num2str(round(space_avg_lambda_list(i), 2)), 'mm'];
            legend(ax, 'interface (OpenFOAM)', 'interface (smoothened)', 'Location', 'south', 'Orientation', 'horizontal');
        end
        
        % removing extra whitespace
        f=getframe(gcf);
        iswhite=min(f.cdata,[],3)==255;
        blankcols=all(iswhite,1);
        col_ind = find(~blankcols,1,'first'):find(~blankcols,1,'last');
        blankrows=all(iswhite,2);
        row_ind = find(~blankrows,1,'first'):find(~blankrows,1,'last');
        cropdata = f.cdata(row_ind,col_ind,:);
        pbaspect([5 1 1]);
        
        %         savename = strcat(saving_path, save_title, '.eps');
%         print(gcf, savename,'-depsc2','-r600');

        savename = strcat(saving_path, save_title, '.png');
%         saveas(ax, savename);
        imwrite(cropdata, savename);
    end
    %---------------------------------------------------

    % finding times correspinding to maximum of peak y coordinates
    % withing 900 and 1000 ms range to be used for avg ambda and attenuation coeff calcluation 
    % or use below filter to get indices of time corresponfing to high amp/ in phase only

    [~, peak_indices_time] = sort(space_max_amp_list);

    % with all peaks for fit and different amp threshold for 1e-1, we can
    % get good results for with and without standing wave cases
    peak_indices_time_sorted = peak_indices_time(end-top_num:end);         
    
    % calculating average parameters
    time_avg_lambda = mean(space_avg_lambda_list(peak_indices_time_sorted), 'omitnan');
    time_avg_k_real = mean(space_avg_k_list(peak_indices_time_sorted), 'omitnan');
    negative_coeff_indices = find(space_attenuation_coeff_list(peak_indices_time_sorted, 2)<0);
    % mean of negative coeffs but used abs value because of e^i(kx-wt) term
    time_avg_k_imag = abs(mean(space_attenuation_coeff_list(peak_indices_time_sorted(negative_coeff_indices), 2), 'omitnan'));

    %---------------------------------------------------
    omega_val = 2*pi*freq;

    sigma_calc = (omega_val^2/time_avg_k_real-g)*rho/time_avg_k_real^2; % N/m, using Eq. 2 from paper, dispersion equation for inviscid liquid
 
    if ~isnan(time_avg_k_imag)
        k_val = time_avg_k_real + 1i*abs(time_avg_k_imag);% simple eq
        mu_calc1 = real(3*rho*(omega_val/real(k_val))*imag(k_val)/4/real(k_val)^2); % using Eq. 3 from paper
        % solving actual dispersion relation (sq root terms) to get single mu values as output
%         mu_calc2 = real(double(vpasolve(sqrt(F(omega_val, k_val, nu, sigma_calc))==0)*rho)); % Pa.s
        mu_calc2 = real(double(vpasolve(F(omega_val, k_val, nu, sigma_calc)==0))*rho); % Pa.s

    else
        mu_calc1 = 0;
        mu_calc2 = 0;
    end
    %---------------------------------------------------
    % writing parameters to file
    
%     filename = strcat(saving_path, '_', name, '_avg_wavelength=', num2str(time_avg_lambda, '%0.1f'), '_mm_avg_k=', num2str(time_avg_k_real, '%0.1f'), '+', num2str(abs(time_avg_k_imag), '%0.1f'), 'i, sigma_c=', num2str(sigma_calc, '%0.1f'), ', mu_c1=', num2str(mu_calc1, '%0.1f'), ', mu_c2=', num2str(mu_calc2, '%0.1f'), '.txt');
    filename = strcat(saving_path, '_', name, '_avg_k=', num2str(time_avg_k_real, '%0.1f'), '+', num2str(abs(time_avg_k_imag), '%0.1f'), 'i, sigma_c=', num2str(sigma_calc*1e3, '%0.1f'), ', mu_c1=', num2str(mu_calc1*1e3, '%0.1f'), ', mu_c2=', num2str(mu_calc2'*1e3, '%0.1f,'), '.txt');

%     if time_avg_k_imag>=0
%         filename = strcat(saving_path, '_', name, '_space-time_avg_lambda=', num2str(time_avg_lambda, '%0.1f'), '_mm_', 'time_avg_complex_k=', num2str(time_avg_k_real, '%0.1f'), '-', num2str(time_avg_k_imag, '%0.1f'), 'i.txt');
%     else
%         filename = strcat(saving_path, '_', name, '_space-time_avg_lambda=', num2str(time_avg_lambda, '%0.1f'), '_mm_', 'time_avg_complex_k=', num2str(time_avg_k_real, '%0.1f'), '+', num2str(abs(time_avg_k_imag), '%0.1f'), 'i.txt');
%     end
    fid = fopen(filename, 'w');
    fprintf(fid, '------------------Top Indices-------------------\n');
    fprintf(fid, 'time\t\tlambda\t\tcomplex wavenumber\n');
    fprintf(fid, '(s)\t\t(mm)\t\t(1/m)\n');
    fprintf(fid, '-----------------------------------------------\n');
    
    for index=1:length(peak_indices_time_sorted)   
        if space_attenuation_coeff_list(peak_indices_time_sorted(index), 2)>=0
            % swapping signs for attenuation coeff before writing in files
            fprintf(fid, '%0.3f\t\t%0.1f\t\t%0.1f-%0.1fi\n', time_steps(peak_indices_time_sorted(index)), space_avg_lambda_list(peak_indices_time_sorted(index)), space_avg_k_list(peak_indices_time_sorted(index)), abs(space_attenuation_coeff_list(peak_indices_time_sorted(index), 2)));
        else
            fprintf(fid, '%0.3f\t\t%0.1f\t\t%0.1f+%0.1fi\n', time_steps(peak_indices_time_sorted(index)), space_avg_lambda_list(peak_indices_time_sorted(index)), space_avg_k_list(peak_indices_time_sorted(index)), abs(space_attenuation_coeff_list(peak_indices_time_sorted(index), 2)));
        end
    end
    fprintf(fid, '-----------------------------------------------\n'); fprintf(fid, '\n'); fprintf(fid, '\n'); fprintf(fid, '\n');
    fprintf(fid, '------------------All Indices------------------\n');
    for index=1:length(time_steps)      
        if space_attenuation_coeff_list(index, 2)>=0
            fprintf(fid, '%0.3f\t\t%0.1f\t\t%0.1f-%0.1fi\n', time_steps(index), space_avg_lambda_list(index), space_avg_k_list(index), abs(space_attenuation_coeff_list(index, 2)));
        else
            fprintf(fid, '%0.3f\t\t%0.1f\t\t%0.1f+%0.1fi\n', time_steps(index), space_avg_lambda_list(index), space_avg_k_list(index), abs(space_attenuation_coeff_list(index, 2)));
        end 
        
        %fprintf(fid, '%0.3f,%0.1f,%0.1f\n', time_steps(index), space_avg_lambda_list(index), space_attenuation_coeff_list(index, 2));
    end    
    fclose('all');

    %---------------------------------------------------
    % plotting 

    figure;
    yyaxis left;
    scatter(time_steps, space_attenuation_coeff_list(:, 2), 'o');
    hold on;
    scatter(time_steps(peak_indices_time_sorted), space_attenuation_coeff_list(peak_indices_time_sorted, 2), 'r*');
    xlabel('time (s)');
    ylabel('exponential fit attenuation coefficient (1/m)');
    grid 'on', grid 'minor';
    yyaxis right;
    plot(time_steps, space_max_amp_list);
    hold on;
    scatter(time_steps(peak_indices_time_sorted), space_max_amp_list(peak_indices_time_sorted), '*');
    xlabel('time (s)');
    ylabel('maximum amplitude (m)');
    grid 'on', grid 'minor';
    savename = strcat(saving_path, '_space_max_amp_list_vs_attenuation_coeff_', name, '.png');
    saveas(gcf, savename); 
    
    figure;
    yyaxis left;
    scatter(time_steps, space_avg_lambda_list(:, 1), 'o');
    hold on;
    scatter(time_steps(peak_indices_time_sorted), space_avg_lambda_list(peak_indices_time_sorted, 1), 'r*');
    xlabel('time (s)');
    ylabel('space averaged wavelength (mm)');
    grid 'on', grid 'minor';
    yyaxis right;
    plot(time_steps, space_max_amp_list);
    hold on;
    scatter(time_steps(peak_indices_time_sorted), space_max_amp_list(peak_indices_time_sorted), '*');
    xlabel('time (s)');
    ylabel('maximum amplitude (m)');
    grid 'on', grid 'minor';
    savename = strcat(saving_path, '_space_max_amp_list_vs_space_avg_lambda_', name, '.png');
    saveas(gcf, savename);  
    close('all');
end
set(0,'DefaultFigureVisible','on');
toc;