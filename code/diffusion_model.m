tic %begin benchmarking

% Diffusion model
% IMPORTANT NOTE: the x-axis and y-axis are reversed due to MATLAB's row-column data storage convention

global x all_cdf diffusivities

%%% PARAMETERS %%%
% Lattice parameters
visualize_lattice = 1; %0: no visualization, 1: visualization
save_lattice = 1; %0: don't save, 1: save (save the lattice to a .mat file if it is newly generated) --> used only if a lattice is generated
heterogeneity = 0.5; %0: perfectly homogeneous, 1: maximal heterogeneity **CURRENTLY NOT USED --> used only if a lattice is generated
lattice_x = 1e4; % --> used only if a lattice is generated
lattice_y = 1e4; % --> used only if a lattice is generated
conversion_factor = 0.1; %conversion factor, units of seconds per time point
single_diffusivity_toggle = 1; %0: multiple subregions, 1: uniform lattice
single_diffusivity = 12500; %value of the diffusivity when constructing a single-diffusivity lattice

% Simulation parameters
time_pts = 5000; %total time points (absolute time, camera frame-rate)
n = 1; %number of simulated particles.
random_start = 0; %0: all particles start at the center of the lattice, 1: particles are each assigned a random start location, -1: all particles start at a hard-coded location

% Plotting parameters
% multiples_delta_time = [1,5,10,50,100,150,200]; %additional time point intervals for displacement histogram generation (each value corresponds to a histogram)
multiples_delta_time = [1,10,100]; %additional time point intervals for displacement histogram generation (each value corresponds to a histogram)

% Workflow parameters
save_data = 1; %0: no data variables are saved, 1: certain data variables are saved (data_matrix, boundary_collision, all_displacement_storage_x, all_displacement_storage_y)

%%% SET-UP %%%
% Fetch a lattice
try
    disp('Attempting to load a pre-existing lattice and associated data.')
    load('lattice.mat','lattice')
    load('lattice_data.mat','x','all_cdf','diffusivities')
    disp('Lattice and associated data loaded from file.')
catch
    disp('No pre-existing lattice is available. Generating a new lattice.')
    tic %begin benchmarking
    [lattice,x,all_cdf,diffusivities] = gen_lattice(save_lattice,heterogeneity,lattice_x,lattice_y,conversion_factor,single_diffusivity_toggle,single_diffusivity);
    toc %end benchmarking
    disp('Lattice generated successfully.')
end

if visualize_lattice == 1
    if length(unique(lattice)) == 1
        lattice_visualization = rescale(lattice,1,1);
    else
        lattice_visualization = rescale(lattice,0,1);
    end
    
    lattice_visualization = fliplr(rot90(lattice_visualization,-1));
    imshow(lattice_visualization)
end

% Failsafe in case a loaded lattice's settings don't match those of the current lattice size
lattice_x = size(lattice,2);
lattice_y = size(lattice,1);

% Matrix to store all relevant simulation data
data_matrix = zeros(n,time_pts,2); %for each particle at each time point, store the current x and y coordinates

%%% WALK SIMULATION %%%
% Memory of boundary collisions
boundary_collision = zeros(1,n);

for i = 1:n %iterate through each particle
    disp_message_particle = strcat(['Simulating particle #' num2str(i) '.']);
    disp(disp_message_particle)
    
    data_matrix(i,1,:) = start_location(random_start,lattice_x,lattice_y); %set the initial position of the particle
    
    for j = 2:time_pts %for each particle, iterate through each time point
        try
            current_diffusivity = lattice(data_matrix(i,j-1,1),data_matrix(i,j-1,2));
        catch
            disp('WARNING. The particle struck the boundary and was rendered immobile.')
            data_matrix(i,j-1,1) = 0; %remove the displacement that crosses the boundary
            data_matrix(i,j-1,2) = 0; %remove the displacement that crosses the boundary
            boundary_collision(i) = 1; %remember that this particle struck the boundary
            break
        end
        
        current_displacement_x = round(get_dispmnt_variation(current_diffusivity)); %randomly select a distance in the x-direction
        current_displacement_y = round(get_dispmnt_variation(current_diffusivity)); %randomly select a distance in the y-direction
        
        % Update the particle's position
        data_matrix(i,j,1) = data_matrix(i,j-1,1) + current_displacement_x;
        data_matrix(i,j,2) = data_matrix(i,j-1,2) + current_displacement_y;

    end
end

if save_data == 1
    save('walk_data.mat','data_matrix','boundary_collision')
end

%%% WALK VISUALIZATION %%%
% Plot the random walk of a single simulated particle
if length(unique(lattice)) == 1
    lattice_visualization = rescale(lattice,1,1);
else
    lattice_visualization = rescale(lattice,0,1);
end

lattice_visualization = fliplr(rot90(lattice_visualization,-1));

for i = 1:n
    figure()
    imshow(lattice_visualization)
    hold on
    
    % Remove excess zeros (in the event the particle struck the boundary)
    current_particle_x = data_matrix(i,:,1);
    current_particle_x(current_particle_x==0) = [];
    
    % Remove excess zeros (in the event the particle struck the boundary)
    current_particle_y = data_matrix(i,:,2);
    current_particle_y(current_particle_y==0) = [];
    
    plot(current_particle_x,current_particle_y) %plot the entire trajectory
    plot(current_particle_x(1),current_particle_y(1),'>g','MarkerFaceColor','g','MarkerSize',10) %mark the start point
    plot(current_particle_x(end),current_particle_y(end),'sr','MarkerFaceColor','r','MarkerSize',10) %mark the end point
    legend('Particle Trajectory','Start','End')
    
    clearvars current_particle_x current_particle_y %clear variables for the next loop iteration
end

%%% HISTOGRAMS %%%
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

% FOR TESTING ("hypotenuse" approach) %

% histogram_data = zeros(n,time_pts,size(multiples_delta_time,2));
% counter_hist = 1;
% for dt = multiples_delta_time
%     for i = 1:n
%         for j = 1:time_pts-dt
%             histogram_data(i,j,counter_hist) = pdist([data_matrix(i,j,1),data_matrix(i,j,2);data_matrix(i,j+dt,1),data_matrix(i,j+dt,2)]);
%             hist_displacement_x = data_matrix(i,j+dt,1) - data_matrix(i,j,1);
%             hist_displacement_y = data_matrix(i,j+dt,2) - data_matrix(i,j,2);
% %             If both displacement components are negative or if the x component is negative and the y component is positive, the displacement is considered to be negative
% %             if (hist_displacement_x<=0 && hist_displacement_y<0) || (hist_displacement_x<=0 && hist_displacement_y>0)
% %                 histogram_data(i,j,counter_hist) = -sqrt(hist_displacement_x^2 + hist_displacement_y^2);
% % %                 histogram_data(i,j,counter_hist) = -(hist_displacement_x^2 + hist_displacement_y^2);
% %             else
% %                 histogram_data(i,j,counter_hist) = sqrt(hist_displacement_x^2 + hist_displacement_y^2);
% % %                 histogram_data(i,j,counter_hist) = (hist_displacement_x^2 + hist_displacement_y^2);
% %             end
% %             histogram_data(i,j,counter_hist) = hist_displacement_y;
% %             histogram_data(i,j,counter_hist) = (hist_displacement_x^2 + hist_displacement_y^2);
%         end
%     end
%     counter_hist = counter_hist + 1;
% end
% 
% % % Remove erroneous displacements that result from a particle striking the boundary
% % for i = 1:n
% %     for j = 1:size(multiples_delta_time,2)
% % %         if boundary_collision(i) == 1 %only modify the displacement data if the given  particle strikes the boundary
% % %             first_zero_idx = find(((histogram_data(i,:,j)==0)+([diff(histogram_data(i,:,j)) 0]==0))==2,1); %find the index of the two consecutive zeros in histogram_data for the given multiple of delta-t
% % %             histogram_data(i,(first_zero_idx-multiples_delta_time(j)):first_zero_idx-1,j) = 0; %remove the appropriate number of erroneous displacements
% % %         end
% %         if boundary_collision(i) == 1 %only modify the displacement data if the given particle strikes the boundary
% %             start_of_trailing_zeros = find(histogram_data(i,:,j),1,'last') + 1;
% %             try
% %                 % Remove the appropriate number of erroneous displacements
% %                 histogram_data(i,(start_of_trailing_zeros-multiples_delta_time(j)):start_of_trailing_zeros-1,j) = 0; 
% %             catch
% %                 % If the particle stopped so early that all displacements are erroneous, set the enter set of displacements to zero
% %                 histogram_data(i,:,j) = 0;
% %             end
% %         end
% %     end
% % end
% 
% % Generate one histogram for each multiple of delta-t
% for i = 1:size(multiples_delta_time,2)
%     clearvars hist_plotting
%     hist_plotting = histogram_data(:,:,i);
%     hist_plotting = reshape(hist_plotting,[1,numel(hist_plotting)]);
%     hist_plotting(hist_plotting==0) = [];
%     hist_plotting = [hist_plotting -hist_plotting];
%     
%     figure()
% %     hist_plotting = rescale(hist_plotting);
% 	hist_obj = histogram(hist_plotting, 'Normalization', 'pdf');
% %     fit_data = fitdist(hist_plotting','Beta'); %obtain fit data
%     fit_data = fitdist(hist_plotting','Normal'); %obtain fit data
% 	eval_vals = (hist_obj.BinEdges(1)-20:0.1:hist_obj.BinEdges(end)+20);
% %     eval_vals = (0:0.001:1);
%     fit_data_pdf = pdf(fit_data,eval_vals); %compute the corresponding PDF
%     hold on
%     plot(eval_vals,fit_data_pdf,'LineWidth',2) %overlay the PDF on top of the histogram
% %     histogram(hist_plotting, 'Normalization', 'probability')
%     title('Step Size Distribution')
%     xlabel('\Deltax, \Deltay [10^{-2}\mum]')
%     hist_y_label_str = strcat(['P(\Deltax, \Deltay, \Delta\tau=' num2str(multiples_delta_time(i)) ' time points)']);
%     ylabel(hist_y_label_str)
%     fit_legend = strcat(['Mean = ' num2str(fit_data.mu) ', Std. Dev. = ' num2str(fit_data.sigma)]);
%     legend('Distribution',fit_legend)
% 
% end

% FOR TESTING ("hypotenuse" approach) %

%%% MSD(DELTA-T) %%%
% MSD(delta-t) plot
squared_displacements_abs_time = zeros(n,time_pts);

% Compute the displacements from the start point to each point in the simulation
for i = 1:n
    for j = 1:time_pts
        displacement_x = data_matrix(i,j,1) - data_matrix(i,1,1);
        displacement_y = data_matrix(i,j,2) - data_matrix(i,1,2);
        squared_displacements_abs_time(i,j) = displacement_x^2 + displacement_y^2;
    end
end

% Generate the squared displacement versus absolute time curves
figure()
for i = 1:n
    plot_sdat = squared_displacements_abs_time(i,:);
    plot_sdat(plot_sdat==plot_sdat(end)) = []; %remove the artifically inflated final values (applicable only when a particle was immobilized by the boundary)
%     plot(plot_sdat,'Color',[0,0,(1/n)*i])
    plot(plot_sdat,'Color',[0,0,0])
    hold on
    
    clearvars plot_sdat
end

% Generate the corresponding MSD(delta-t) curve
msd_dt_data = zeros(size(squared_displacements_abs_time));
for i = 1:n
    plot_msd_dt = squared_displacements_abs_time(i,:);
    plot_msd_dt(plot_msd_dt==plot_msd_dt(end)) = 0; %replace the artifically inflated final values with zero (applicable only when a particle was immobilized by the boundary)
    msd_dt_data(i,:) = plot_msd_dt;

    clearvars plot_msd_dt
end

% Plot the corresponding MSD(delta-t) curve
msd_dt = sum(msd_dt_data,1)./sum(msd_dt_data~=0,1);
plot(msd_dt,'LineWidth',2,'Color','red')
title('Squared Displacement vs. Absolute Time')
xlabel('Absolute Time \Deltat [simulation time points]')
ylabel('MSD(\Deltat) [(10^{-2}\mum)^2]')

%%% TIME-AVERAGED MSD(DELTA-TAU) %%%
% Time-averaged MSD(delta-tau) plot
delta_taus = 1:(time_pts/10); %time point intervals for displacement measurements (given by the first 10% of time points)
sqd_dispmnts_lag_time = zeros(n,time_pts,size(delta_taus,2)); %storage for displacements at each time point interval

% Compute the displacements for the given delta-tau values (multiples of delta-t)
counter_msd_tau = 1;
for dt = delta_taus
    for i = 1:n
        for j = 1:time_pts-dt
            tamsd_plot_displacement_x = data_matrix(i,j+dt,1) - data_matrix(i,j,1);
            tamsd_plot_displacement_y = data_matrix(i,j+dt,2) - data_matrix(i,j,2);
            sqd_dispmnts_lag_time(i,j,counter_msd_tau) = tamsd_plot_displacement_x^2 + tamsd_plot_displacement_y^2;
        end
    end
    counter_msd_tau = counter_msd_tau + 1;
end

% Remove erroneous displacements that result from a particle striking the boundary
for i = 1:n
    for j = 1:size(delta_taus,2)
        if boundary_collision(i) == 1 %only modify the displacement data if the given particle strikes the boundary
            start_of_trailing_zeros = find(sqd_dispmnts_lag_time(i,:,j),1,'last') + 1;
            try
                % Remove the appropriate number of erroneous displacements
                sqd_dispmnts_lag_time(i,(start_of_trailing_zeros-delta_taus(j)):start_of_trailing_zeros-1,j) = 0; 
            catch
                % If the particle stopped so early that all displacements are erroneous, set the enter set of displacements to zero
                sqd_dispmnts_lag_time(i,:,j) = 0;
            end
        end
    end
end

% Set up the plot for the time-averaged squared displacements versus lag time curves
sdlt_plotting = zeros(n,size(delta_taus,2));
for i = 1:n
    for j = 1:size(delta_taus,2)
        current_particle_tau = sqd_dispmnts_lag_time(i,:,j); %isolate one particle
        start_of_trailing_zeros = find(sqd_dispmnts_lag_time(i,:,j),1,'last') + 1; %find the index of the two consecutive zeros in sqd_dispmnts_lag_time for the given multiple of delta-t
        current_particle_tau(start_of_trailing_zeros:end) = []; %remove all trailing zeros from the data matrix
        
        sdlt_plotting(i,j) = mean(current_particle_tau); %compute the mean for the given particle at the given time lag
    end
end

% Plot the time-averaged squared displacements versus lag time curves
figure()
for i = 1:n
    plot(sdlt_plotting(i,:))
    hold on
end

% Convert NaN values (resulting from entirely-zero walks) to zeros to avoid breaking the MSD(delta-tau) curve
sdlt_plotting(isnan(sdlt_plotting)) = 0;

% Set up the plot for the corresponding time-averaged MSD(delta-tau) curve
msd_tau_plotting = zeros(1,size(delta_taus,2));
for i = 1:size(delta_taus,2)
    sdlt_plotting_no_zeros = sdlt_plotting(:,i);
    sdlt_plotting_no_zeros(sdlt_plotting_no_zeros==0) = [];
    msd_tau_plotting(i) = mean(sdlt_plotting_no_zeros);
end

% Plot the corresponding time-averaged MSD(delta-tau) curve
plot(msd_tau_plotting,'LineWidth',2,'Color','red')
title('Time-Averaged Squared Displacement vs. Lag Time')
xlabel('Lag Time \Delta\tau [simulation time points]')
ylabel('Time-Averaged MSD(\Delta\tau) [(10^{-2}\mum)^2]')

% Fetch the line of best fit for the MSD(tau)
msd_tau_fit = polyfit(log10(delta_taus(1:end/2)),log10(msd_tau_plotting(1:end/2)),1);
disp(msd_tau_fit)

toc %end benchmarking