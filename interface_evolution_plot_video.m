clear; clc;
close all;

folder_path = ['C:\Users\mukim\Documents\PhD\Work\Projects\Capillary Wave Surface Reconstruction\Journal Paper 3\'];
parameter = ['pl_vs_ax_dataset\pl\sigma\0.50_20'];
saving_path = strcat(folder_path, parameter, '\processed_data\');
if ~isfolder(saving_path)
    disp('processed folder not found...creating the folder');
    mkdir(saving_path);
end


% write plot in video file
par_index = strfind(parameter, '\');
name = parameter(par_index(1)+1:end);
name = strrep(name,'\','_');

if isfile(strcat(saving_path, name, '.avi'))
    % delete existing file before running it again
   delete(strcat(saving_path, name, '.avi'));
end

axisymm_flag = length(strfind(parameter, '\axisymm\'));

deltaT = 0.001;

if ~contains(parameter, '\mu\')
    time_steps = 950:1:1000;
else
    time_steps = 950:1:1000;
end

filepath = strcat(saving_path, name, '.avi');
if isfile(filepath)
    disp('file found...deleting it before running the code again');
    delete(filepath);
end
vid = VideoWriter(filepath);
vid.FrameRate = 1;
open(vid);

%-------------------------------------------------------------
% plot
fig = figure;
fig.Position = [30 30 1200 400];
ax = axes(fig);
% common_x_axis = 5:1e-2:75; % common x axis for overlapped interface
% common_y_axis = zeros(length(common_x_axis), 1);
% ft = fittype('smoothingspline');
% opts = fitoptions('Method', 'SmoothingSpline');
% opts.SmoothingParam = 0.99;  


for i=1:length(time_steps)    
    time = time_steps(i)*deltaT;
%     disp(time);
%     cla(ax);
    time_file = strcat(folder_path, parameter, '\', num2str(time), '\freeSurf\points');
%     disp(time_file);
    fid = fopen(time_file, 'r');
    tline = fgetl(fid);
    data = [];
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
    
    x = data(1, :)*1e3;
    y = data(2, :)*1e3-10.0;
    scatter(ax, x, y, '.');
    hold (ax, 'on');
    if axisymm_flag
        xlabel(ax, 'radius (mm)');
    else
        xlabel(ax, 'length (mm)');   
    end
    ylabel(ax, 'amplitude (mm)');
    % pbaspect([6 1 1]);

    xlim(ax,[0 80]);   
    if axisymm_flag
        ylim(ax,[-0.05 0.05]);   
    else
        ylim(ax,[-0.05 0.05]);   

%         ylim(ax,[-0.2 0.2]); 
    end
    grid (ax, 'on');
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';

    plot_title = ['time=', num2str(time*1000), 'ms'];
    title(plot_title);  
    set(ax,'FontSize', 10,'FontName', 'Times');

    F = getframe(fig);
    writeVideo(vid,F);
%     pause(0.5);
    
%     % extract overlapped interface
%     
%     [fitresult, ~] = fit(x', y', ft, opts);
%     fitted_y_axis = fitresult(common_x_axis);
%     indices = find(fitted_y_axis > common_y_axis);
%     common_y_axis(indices) = fitted_y_axis(indices);

end
close(vid);

% [fitresult, ~] = fit(common_x_axis', common_y_axis, ft, opts); 
% common_y_axis = fitresult(common_x_axis); % smoothening
% 
% plot(common_x_axis, common_y_axis, 'b');
% hold on;
% % [y_peak_list, x_peak_list] = findpeaks(common_y_axis, common_x_axis, 'MinPeakHeight', 0.1*max(common_y_axis));
% 
% 
% % First derivative
% dy = diff(common_y_axis);
% % Second derivative
% d2y = diff(common_y_axis, 2);
% % Find zero-crossings in dy (from + to -) and d2y < 0 for peaks
% valid_peaks = find(dy(1:end-1).*dy(2:end)<=0 & d2y(1:end)<0);
% 
% peak_indices = valid_peaks + 1;
% valid_peak_indices = peak_indices(find(common_y_axis(peak_indices)>0.1*max(common_y_axis)));
% x_peak_list = common_x_axis(valid_peak_indices);
% y_peak_list = common_y_axis(valid_peak_indices);
% 
% scatter(x_peak_list, y_peak_list, 'b*');
% hold on;
% [fitresult, ~] = fit(x_peak_list', y_peak_list, 'exp1'); 
% fit_common_y_axis = fitresult(common_x_axis);
% plot(common_x_axis, fit_common_y_axis, 'r--');
% 
% disp(round(fitresult.b*-1*1000, 2));