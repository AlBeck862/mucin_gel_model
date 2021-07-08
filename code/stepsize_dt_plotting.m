function stepsize_dt_plotting(n,time_pts,data_matrix,lattice,moving_avg_kernel,conversion_factor,multiplier,diffusivities)
% STEPSIZE_IN_TIME Plot the step size over time for each particle in both
% the x-direction and the y-direction. Plot the frequency of the
% diffusivities encountered by each particle throughout the simulation.

% Time in seconds used for more realistic plotting
rescaled_time = conversion_factor:conversion_factor:time_pts*conversion_factor;

% Storage for the diffusivities encountered by all the particles
all_particle_diffs = cell(n,1);

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
    
    current_particle_diffs = current_particle_diffs(1:end-1);   %match current_particle_diffs with stepsize_x and stepsize_y (the final diffusivity does not matter)
	[unique_diffs,~,idx_c] = unique(current_particle_diffs);    %determine the unique diffusivities encountered
    tally = accumarray(idx_c,1);                                %count the number of times each unique diffusivity occurs
    tally = tally';
    
    all_particle_diffs{i} = current_particle_diffs;           %collect the data from each particle for the heatmap
    
    % Bar plot of the diffusivities encountered by the current particle
    figure()
    unique_diffs = unique_diffs./multiplier;
    bar_data = bar(unique_diffs,tally,0.5);
    xtips = bar_data(1).XEndPoints;
    ytips = bar_data(1).YEndPoints;
    labels = string(bar_data(1).YData);
    text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
    diff_title_str = strcat(['Frequency of Diffusivities (Particle #' num2str(i) ')']);
    title(diff_title_str)
    xlabel('Diffusivity [\mum^2/s]')
    ylabel('Frequency')
    grid on
    
	file_str = strcat(['/temp_results/diffusivities/diffusivities_particle' num2str(i) '.jpeg']);
    saveas(gcf,[pwd file_str]);
    
    % Plot of the x-direction step-sizes of the current particle (corresponding diffusivities are overlaid)
    figure()
    plot(rescaled_time(1:length(stepsizes_x)),movmean(stepsizes_x./(sqrt(multiplier)),moving_avg_kernel)) %plot using a moving average to smooth out the curve
    x_dir_title_str = strcat(['Step-Size (x-direction) and Diffusivity vs. Absolute Time (Particle #' num2str(i) ')']);
    title(x_dir_title_str)
    xlabel('Absolute Time [seconds]')
    ylabel('Step-Size |x| [\mum]')
    
    yyaxis right
    plot(rescaled_time(1:length(current_particle_diffs)),current_particle_diffs./multiplier)
%     plot(movmean(current_particle_diffs,moving_avg_kernel)) %plot using a moving average to smooth out the curve
    ylabel('Diffusivity [\mum^2/s]')
    
	file_str = strcat(['/temp_results/step_sizes/stepsizes_dx_particle' num2str(i) '.jpeg']);
    saveas(gcf,[pwd file_str]);
    
    % Plot of the y-direction step-sizes of the current particle (corresponding diffusivities are overlaid)
    figure()
    plot(rescaled_time(1:length(stepsizes_y)),movmean(stepsizes_y./(sqrt(multiplier)),moving_avg_kernel)) %plot using a moving average to smooth out the curve
    y_dir_title_str = strcat(['Step-Size (y-direction) and Diffusivity vs. Absolute Time (Particle #' num2str(i) ')']);
    title(y_dir_title_str)
    xlabel('Absolute Time [seconds]')
    ylabel('Step-Size |y| [\mum]')
    
	yyaxis right
    plot(rescaled_time(1:length(current_particle_diffs)),current_particle_diffs./multiplier)
%     plot(movmean(current_particle_diffs,moving_avg_kernel)) %plot using a moving average to smooth out the curve
    ylabel('Diffusivity [\mum^2/s]')
    
	file_str = strcat(['/temp_results/step_sizes/stepsizes_dy_particle' num2str(i) '.jpeg']);
    saveas(gcf,[pwd file_str]);
    
end

% Heatmap of the diffusivities encountered by every particle
possible_diffusivities = [diffusivities;diffusivities(end)+1]; %append an extra value for compatibility with the groupcounts function
all_particle_diffs_counts = zeros(i,length(possible_diffusivities)-1);
for i = 1:i
    [single_particle_diffs_counts,~] = groupcounts(all_particle_diffs{i}',possible_diffusivities','IncludeEmptyGroups',true);
    all_particle_diffs_counts(i,:) = single_particle_diffs_counts;
end

h = heatmap(all_particle_diffs_counts);
% h.XDisplayLabels = diffs;
% h.YDisplayLabels = particles;

file_str = strcat('/temp_results/diffusivity_heatmap.jpeg');
saveas(gcf,[pwd file_str]);

end