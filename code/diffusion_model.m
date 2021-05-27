tic

% Diffusion model
% IMPORTANT NOTE: the x- and y-axis are reversed due to MATLAB row-column data storage convention

global x all_cdf diffusivities

%%% PARAMETERS %%%
% Lattice parameters
visualize_lattice = 1; %0: no visualization, 1: visualization
save_lattice = 1; %0: don't save, 1: save (save the lattice to a .mat file if it is newly generated) --> used only if a lattice is generated
heterogeneity = 0.5; %0: perfectly homogeneous, 1: maximal heterogeneity **CURRENTLY NOT USED --> used only if a lattice is generated
lattice_x = 1e4; % --> used only if a lattice is generated
lattice_y = 1e4; % --> used only if a lattice is generated
tau = 0.1; %step unit (this must still be related back to time) --> used only if a lattice is generated

% Simulation parameters
time_pts = 5000; %total time points (absolute time, camera frame-rate)
n = 25; %number of simulated particles.
random_start = 1; %0: all particles start at the center of the lattice, 1: particles are each assigned a random start location

% Plotting parameters
multiples_delta_time = [5,10,50,100]; %additional time point intervals for displacement histogram generation (each value corresponds to a histogram)

%%% SET UP %%%
% Fetch a lattice
try
    disp('Attempting to load a pre-existing lattice and associated data.')
    load('lattice.mat','lattice')
    load('lattice_data.mat','x','all_cdf','diffusivities')
    disp('Lattice and associated data loaded from file.')
catch
    disp('No pre-existing lattice is available. Generating a new lattice.')
    tic %begin benchmarking
    [lattice,x,all_cdf,diffusivities] = gen_lattice(save_lattice,heterogeneity,lattice_x,lattice_y,tau);
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

%%% PARTICLE START LOCATIONS %%%
% Particle start location
if random_start == 1
    % A margin must be set to avoid starting too close to the boundary
    margin_percent_x = 10; %the particle won't start within 10% of the lattice's x-length from the boundary on both sides (e.g.: lattice_x=10000 --> the particle's starting x-coordinate can be between 1000 and 9000)
    start_margin_x = round(lattice_x/margin_percent_x);
    min_start_x = start_margin_x;
    max_start_x = lattice_x-start_margin_x;
    
    % A margin must be set to avoid starting too close to the boundary
    margin_percent_y = 10; %the particle won't start within 10% of the lattice's y-length from the boundary on both sides (e.g.: lattice_y=10000 --> the particle's starting y-coordinate can be between 1000 and 9000)
    start_margin_y = round(lattice_y/margin_percent_y);
    min_start_y = start_margin_y;
    max_start_y = lattice_y-start_margin_y;
    
    % Generate a random start location for each particle
    for i = 1:n
        % Select a random start location at a sufficient distance from the boundary
        x_start = round(min_start_x + (max_start_x-min_start_x).*rand());
        y_start = round(min_start_y + (max_start_y-min_start_y).*rand());
        
        % Set up particle start locations and initial diffusivities
        data_matrix(i,1,1) = x_start;
        data_matrix(i,1,2) = y_start;
    end
elseif random_start == -1 %for testing purposes
    %  The particle will start at a specific, hard-coded location
	x_start = 200;
    y_start = 200;
    
    % Set up particle start locations and initial diffusivities
	data_matrix(:,1,1) = x_start;
    data_matrix(:,1,2) = y_start;
else
    % The particle will start at the center of the lattice
    x_start = round(lattice_x/2);
    y_start = round(lattice_y/2);
    
    % Set up particle start locations and initial diffusivities
    data_matrix(:,1,1) = x_start;
    data_matrix(:,1,2) = y_start;
end

%%% WALK SIMULATION %%%
% Storage for the global histogram of displacements
current_displacement_storage = zeros(n,time_pts-1);

% Memory of boundary collisions
boundary_collision = zeros(1,n);

for i = 1:n %iterate through each particle
    disp_message_particle = strcat(['Simulating particle #' num2str(i) '.']);
    disp(disp_message_particle)
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
        
        current_displacement = round(get_dispmnt_variation(current_diffusivity)); %randomly select a distance
        
        direction_select = randi(4); %randomly select a direction
        
        if direction_select == 1 %+x (RIGHT)
            data_matrix(i,j,1) = data_matrix(i,j-1,1) + abs(current_displacement);
            data_matrix(i,j,2) = data_matrix(i,j-1,2);
            current_displacement_storage(i,j-1) = current_displacement;
        elseif direction_select == 2 %-x (LEFT)
            data_matrix(i,j,1) = data_matrix(i,j-1,1) - abs(current_displacement);
            data_matrix(i,j,2) = data_matrix(i,j-1,2);
            current_displacement_storage(i,j-1) = -current_displacement;
        elseif direction_select == 3 %+y (DOWN)
            data_matrix(i,j,1) = data_matrix(i,j-1,1);
            data_matrix(i,j,2) = data_matrix(i,j-1,2) + abs(current_displacement);
            current_displacement_storage(i,j-1) = -current_displacement;
        elseif direction_select == 4 %-y (UP)
            data_matrix(i,j,1) = data_matrix(i,j-1,1);
            data_matrix(i,j,2) = data_matrix(i,j-1,2) - abs(current_displacement);
            current_displacement_storage(i,j-1) = current_displacement;
        end
    end
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
% Histogram for all displacements (1*(delta-t))
clearvars all_displacement_storage
all_displacement_storage = [];
for i = 1:n
    % Store the processed data for histogram plotting
    all_displacement_storage = [all_displacement_storage no_trailing_zeros(current_displacement_storage(i,:))];
end

figure()
histogram(all_displacement_storage, 'Normalization', 'pdf')
title('Step Size Distribution')
xlabel('\Deltax, \Deltay [10^{-2}\mum]')
ylabel('P(\Deltax, \Deltay, \Delta\tau=1 time point)')

% Histograms for displacements using greater multiples of delta-t (iteratively: multiples_delta_time*(delta-t))
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
    clearvars all_displacement_storage_x all_displacement_storage_y all_displacement_storage_both
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
        
        % FIRST METHOD (FILTERING) %
        current_histogram_data_x = histogram_data_x(i,:,j);
        current_histogram_data_y = histogram_data_y(i,:,j);
        
        current_histogram_data_x_copy = histogram_data_x(i,:,j);
        current_histogram_data_y_copy = histogram_data_y(i,:,j);

        current_histogram_data_x((current_histogram_data_x_copy==0 & current_histogram_data_y_copy~=0)) = [];
        current_histogram_data_y((current_histogram_data_y_copy==0 & current_histogram_data_x_copy~=0)) = [];
        
        current_histogram_data_x = no_trailing_zeros(current_histogram_data_x);
        current_histogram_data_y = no_trailing_zeros(current_histogram_data_y);
        
        % Store the processed data for histogram plotting
        all_displacement_storage_x = [all_displacement_storage_x current_histogram_data_x];
        all_displacement_storage_y = [all_displacement_storage_y current_histogram_data_y];
        % FIRST METHOD (FILTERING) %
        
        % SECOND METHOD (NO FILTERING) %
        % Store the processed data for histogram plotting
%         all_displacement_storage_x = [all_displacement_storage_x no_trailing_zeros(histogram_data_x(i,:,j))];
%         all_displacement_storage_y = [all_displacement_storage_y no_trailing_zeros(histogram_data_y(i,:,j))];
        % SECOND METHOD (NO FILTERING) %
        
    end
    
    % Combine the direction-specific data for histogram plotting
    all_displacement_storage_both = [all_displacement_storage_x all_displacement_storage_y];
    
%     figure()
%     histogram(all_displacement_storage_x, 'Normalization', 'pdf')
%     title('Step Size Distribution')
%     xlabel('\Deltax [10^{-2}\mum]')
%     hist_y_label_str = strcat(['P(\Deltax, \Delta\tau=' num2str(multiples_delta_time(j)) ' time points)']);
%     ylabel(hist_y_label_str)
%     
%     figure()
%     histogram(all_displacement_storage_y, 'Normalization', 'pdf')
%     title('Step Size Distribution')
%     xlabel('\Deltay [10^{-2}\mum]')
%     hist_y_label_str = strcat(['P(\Deltay, \Delta\tau=' num2str(multiples_delta_time(j)) ' time points)']);
%     ylabel(hist_y_label_str)

	figure()
    histogram(all_displacement_storage_both, 'Normalization', 'pdf')
    title('Step Size Distribution')
    xlabel('\Deltax, \Deltay [10^{-2}\mum]')
    hist_y_label_str = strcat(['P(\Deltax, \Deltay, \Delta\tau=' num2str(multiples_delta_time(j)) ' time points)']);
    ylabel(hist_y_label_str)

end

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
        if boundary_collision(i) == 1 %only modify the displacement data if the given  particle strikes the boundary
            first_zero_idx = find(((sqd_dispmnts_lag_time(i,:,j)==0)+([diff(sqd_dispmnts_lag_time(i,:,j)) 0]==0))==2,1); %find the index of the two consecutive zeros in sqd_dispmnts_lag_time for the given multiple of delta-t
            try
                sqd_dispmnts_lag_time(i,(first_zero_idx-delta_taus(j)):first_zero_idx-1,j) = 0; %remove the appropriate number of erroneous displacements
            catch
                sqd_dispmnts_lag_time(i,:,j) = 0; %remove all displacements if the particle is immobilized at a time point lesser than the value of the time lag
            end
        end
    end
end

% Set up the plot for the time-averaged squared displacements versus lag time curves
sdlt_plotting = zeros(n,size(delta_taus,2));
for i = 1:n
    for j = 1:size(delta_taus,2)
        current_particle_tau = sqd_dispmnts_lag_time(i,:,j); %isolate one particle
        first_zero_idx = find(((sqd_dispmnts_lag_time(i,:,j)==0)+([diff(sqd_dispmnts_lag_time(i,:,j)) 0]==0))==2,1); %find the index of the two consecutive zeros in sqd_dispmnts_lag_time for the given multiple of delta-t
        current_particle_tau(first_zero_idx:end) = []; %remove all trailing zeros from the data matrix
        
        sdlt_plotting(i,j) = mean(current_particle_tau); %compute the mean for the given particle at the given time lag
    end
end

% Plot the time-averaged squared displacements versus lag time curves
figure()
for i = 1:n
    plot(sdlt_plotting(i,:))
    hold on
end

% FIRST METHOD %
% % Remove trajectories that would break the MSD(delta-tau) curve
% while any(any(isnan(sdlt_plotting))) %loop as long as there are NaN values in the matrix
%     for i = 1:size(sdlt_plotting,1) %loop over the current size of the matrix (since rows may be removed, the matrix's size may change with each while loop iteration)
%         if any(isnan(sdlt_plotting(i,:)))
%             sdlt_plotting(i,:) = []; %remove rows (particle trajectories) that contain NaN values
%             break
%         end
%     end
% end
% 
% % Set up the plot for the corresponding time-averaged MSD(delta-tau) curve
% msd_tau_plotting = mean(sdlt_plotting,1);
% FIRST METHOD

% SECOND METHOD %
% Convert NaN values (resulting from entirely-zero walks) to zeros to avoid breaking the MSD(delta-tau) curve
sdlt_plotting(isnan(sdlt_plotting)) = 0;

% Set up the plot for the corresponding time-averaged MSD(delta-tau) curve
msd_tau_plotting = zeros(1,size(delta_taus,2));
for i = 1:size(delta_taus,2)
    sdlt_plotting_no_zeros = sdlt_plotting(:,i);
    sdlt_plotting_no_zeros(sdlt_plotting_no_zeros==0) = [];
    msd_tau_plotting(i) = mean(sdlt_plotting_no_zeros);
end
% SECOND METHOD %

% Plot the corresponding time-averaged MSD(delta-tau) curve
plot(msd_tau_plotting,'LineWidth',2,'Color','red')
title('Time-Averaged Squared Displacement vs. Lag Time')
xlabel('Lag Time \Delta\tau [simulation time points]')
ylabel('Time-Averaged MSD(\Delta\tau) [(10^{-2}\mum)^2]')

toc