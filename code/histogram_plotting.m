function histogram_plotting(n,time_pts,multiples_delta_time,data_matrix,boundary_collision,save_data)
% HISTOGRAM_PLOTTING Plots two histograms (x-direction displacements and
% y-direction displacements) for each multiple of delta-time.

% Histograms for displacements using predefined multiples of delta-t (iteratively: multiples_delta_time*(delta-t))
histogram_data_x = zeros(n,time_pts,size(multiples_delta_time,2));
histogram_data_y = zeros(n,time_pts,size(multiples_delta_time,2));
counter_hist = 1;
for dt = multiples_delta_time
    for i = 1:n
        for j = 1:time_pts-dt
            histogram_data_x(i,j,counter_hist) = data_matrix(i,j+dt,1) - data_matrix(i,j,1);
            histogram_data_y(i,j,counter_hist) = data_matrix(i,j+dt,2) - data_matrix(i,j,2);
        end
    end
    counter_hist = counter_hist + 1;
end

for j = 1:size(multiples_delta_time,2)
    clearvars all_displacement_storage_x all_displacement_storage_y all_displacement_storage_both eval_vals hist_obj fit_data_pdf
    all_displacement_storage_x = [];
    all_displacement_storage_y = [];
    for i = 1:n
        % Remove erroneous displacements that result from a particle striking the boundary
        if boundary_collision(i) == 1 %only modify the displacement data if the given particle strikes the boundary
            start_of_trailing_zeros = find(histogram_data_x(i,:,j),1,'last') + 1;
            try
                % Remove the appropriate number of erroneous displacements
                histogram_data_x(i,(start_of_trailing_zeros-multiples_delta_time(j)):start_of_trailing_zeros-1,j) = 0; 
                histogram_data_y(i,(start_of_trailing_zeros-multiples_delta_time(j)):start_of_trailing_zeros-1,j) = 0;
            catch
                % If the particle stopped so early that all displacements are erroneous, set the enter set of displacements to zero
                histogram_data_x(i,:,j) = 0;
                histogram_data_y(i,:,j) = 0;
            end
        end
        
        % Store the processed data for histogram plotting
        all_displacement_storage_x = [all_displacement_storage_x no_trailing_zeros(histogram_data_x(i,:,j))];
        all_displacement_storage_y = [all_displacement_storage_y no_trailing_zeros(histogram_data_y(i,:,j))];
        
    end
    
    if save_data == 1
        save_name_str = strcat(['histograms_all_displacements_' num2str(multiples_delta_time(j)) 'dt.mat']);
        save(save_name_str,'all_displacement_storage_x','all_displacement_storage_y')
    end
    
    % Combine the direction-specific data for histogram plotting
    all_displacement_storage_both = [all_displacement_storage_x all_displacement_storage_y];
    
    % x-direction displacements histogram
    figure()
    hist_obj_x = histogram(all_displacement_storage_x, 'Normalization', 'pdf');
    fit_data_x = fitdist(all_displacement_storage_x','Normal'); %obtain fit data
    
	eval_vals_x = (hist_obj_x.BinEdges(1)-20:0.1:hist_obj_x.BinEdges(end)+20);
    fit_data_pdf_x = pdf(fit_data_x,eval_vals_x); %compute the corresponding PDF
    hold on
    plot(eval_vals_x,fit_data_pdf_x,'LineWidth',2) %overlay the PDF on top of the histogram
    
    title('Step Size Distribution')
    xlabel('\Deltax [10^{-2}\mum]')
    hist_y_label_str = strcat(['P(\Deltax, \Delta\tau=' num2str(multiples_delta_time(j)) ' time points)']);
    ylabel(hist_y_label_str)
	fit_legend_x = strcat(['Mean = ' num2str(fit_data_x.mu) ', Std. Dev. = ' num2str(fit_data_x.sigma)]);
    legend('Distribution',fit_legend_x)
    
	file_str = strcat(['/temp_results/histograms/histogram_' num2str(multiples_delta_time(j)) 'dt_dx.jpeg']);
    saveas(gcf,[pwd file_str]);
    
	% y-direction displacements histogram
    figure()
    hist_obj_y = histogram(all_displacement_storage_y, 'Normalization', 'pdf');
	fit_data_y = fitdist(all_displacement_storage_y','Normal'); %obtain fit data
    
	eval_vals_y = (hist_obj_y.BinEdges(1)-20:0.1:hist_obj_y.BinEdges(end)+20);
    fit_data_pdf_y = pdf(fit_data_y,eval_vals_y); %compute the corresponding PDF
    hold on
    plot(eval_vals_y,fit_data_pdf_y,'LineWidth',2) %overlay the PDF on top of the histogram
    
    title('Step Size Distribution')
    xlabel('\Deltay [10^{-2}\mum]')
    hist_y_label_str = strcat(['P(\Deltay, \Delta\tau=' num2str(multiples_delta_time(j)) ' time points)']);
    ylabel(hist_y_label_str)
	fit_legend_y = strcat(['Mean = ' num2str(fit_data_y.mu) ', Std. Dev. = ' num2str(fit_data_y.sigma)]);
    legend('Distribution',fit_legend_y)

	file_str = strcat(['/temp_results/histograms/histogram_' num2str(multiples_delta_time(j)) 'dt_dy.jpeg']);
    saveas(gcf,[pwd file_str]);
    
    % Stacked displacements histogram
% 	figure()
%     hist_obj = histogram(all_displacement_storage_both, 'Normalization', 'pdf');
%     fit_data = fitdist(all_displacement_storage_both','Normal'); %obtain fit data
%     
% %     fit_data = fitdist(all_displacement_storage_both','Logistic');
% %     fit_data = fitdist(all_displacement_storage_both','Stable');
%     
%     eval_vals = (hist_obj.BinEdges(1)-20:0.1:hist_obj.BinEdges(end)+20);
%     fit_data_pdf = pdf(fit_data,eval_vals); %compute the corresponding PDF
%     hold on
%     plot(eval_vals,fit_data_pdf,'LineWidth',2) %overlay the PDF on top of the histogram
%     
%     title('Step Size Distribution')
%     xlabel('\Deltax, \Deltay [10^{-2}\mum]')
%     hist_y_label_str = strcat(['P(\Deltax, \Deltay, \Delta\tau=' num2str(multiples_delta_time(j)) ' time points)']);
%     ylabel(hist_y_label_str)
%     fit_legend = strcat(['Mean = ' num2str(fit_data.mu) ', Std. Dev. = ' num2str(fit_data.sigma)]);
%     legend('Distribution',fit_legend)
    
end

end