function stepsize_dt_plotting(n,data_matrix,lattice,moving_avg_kernel)
% STEPSIZE_IN_TIME Plot the step size over time for each particle in both
% the x-direction and the y-direction. Plot the frequency of the
% diffusivities encountered by each particle throughout the simulation.

for i = 1:n
    % Isolate the trajectory of the current particle
    current_particle_x = no_trailing_zeros(data_matrix(i,:,1));
    current_particle_y = no_trailing_zeros(data_matrix(i,:,2));
    
    % Compute the distances between each position
    stepsizes_x = abs(diff(current_particle_x));
    stepsizes_y = abs(diff(current_particle_y));
    
    % Set the total number of time points to be the length of the shortest data array
    if length(current_particle_x) == length(current_particle_y)
        time_pts = length(current_particle_x);
    elseif length(current_particle_x) > length(current_particle_y)
        time_pts = length(current_particle_y);
    else
        time_pts = length(current_particle_x);
    end
        
    % Collect the diffusivities the current particle encountered
    current_particle_diffs = zeros(1,time_pts);
    for j = 1:time_pts
        current_particle_diffs(j) = lattice(current_particle_x(j),current_particle_y(j));
    end
    
    current_particle_diffs = current_particle_diffs(1:end-1); %match current_particle_diffs with stepsize_x and stepsize_y (the final diffusivity does not matter)
	[unique_diffs,~,idx_c] = unique(current_particle_diffs); %determine the unique diffusivities encountered
    tally = accumarray(idx_c,1); %count the number of times each unique diffusivity occurs
    tally = tally';
    
    % Bar plot of the diffusivities encountered by the current particle
    figure()
%     histogram(current_particle_diffs)
    bar_data = bar(unique_diffs,tally,0.5);
    xtips = bar_data(1).XEndPoints;
    ytips = bar_data(1).YEndPoints;
    labels = string(bar_data(1).YData);
    text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
    diff_title_str = strcat(['Frequency of Diffusivities (Particle #' num2str(i) ')']);
    title(diff_title_str)
    xlabel('Diffusivity [(10^{-4}\mum)^2/s]')
    ylabel('Frequency')
    grid on
    
	file_str = strcat(['/temp_results/diffusivities/diffusivities_particle' num2str(i) '.jpeg']);
    saveas(gcf,[pwd file_str]);
    
    % Plot of the x-direction step-sizes of the current particle (corresponding diffusivities are overlaid)
    figure()
    plot(movmean(stepsizes_x,moving_avg_kernel))
    x_dir_title_str = strcat(['Step-Size (x-direction) and Diffusivity vs. Absolute Time (Particle #' num2str(i) ')']);
    title(x_dir_title_str)
    xlabel('Absolute Time [time points]')
    ylabel('Step-Size |x| [(10^{-2}\mum)]')
    
    yyaxis right
    plot(current_particle_diffs)
%     plot(movmean(current_particle_diffs,moving_avg_kernel))
    ylabel('Diffusivity [(10^{-4}\mum)^2/s]')
    
	file_str = strcat(['/temp_results/step_sizes/stepsizes_dx_particle' num2str(i) '.jpeg']);
    saveas(gcf,[pwd file_str]);
    
    % Plot of the y-direction step-sizes of the current particle (corresponding diffusivities are overlaid)
    figure()
    plot(movmean(stepsizes_y,moving_avg_kernel))
    y_dir_title_str = strcat(['Step-Size (y-direction) and Diffusivity vs. Absolute Time (Particle #' num2str(i) ')']);
    title(y_dir_title_str)
    xlabel('Absolute Time [time points]')
    ylabel('Step-Size |y| [(10^{-2}\mum)]')
    
	yyaxis right
    plot(current_particle_diffs)
%     plot(movmean(current_particle_diffs,moving_avg_kernel))
    ylabel('Diffusivity [(10^{-4}\mum)^2/s]')
    
	file_str = strcat(['/temp_results/step_sizes/stepsizes_dy_particle' num2str(i) '.jpeg']);
    saveas(gcf,[pwd file_str]);
    
end

end