% Diffusion model
% IMPORTANT NOTE: the x- and y-axis are reversed due to MATLAB row-column data storage convention

%%% PARAMETERS %%%
% Lattice parameters
visualize_lattice = 1; %0: no visualization, 1: visualization
save_lattice = 1; %0: don't save, 1: save (save the lattice to a .mat file if it is newly generated)
heterogeneity = 0.5; %0: perfectly homogeneous, 1: maximal heterogeneity **CURRENTLY NOT USED
lattice_x = 1e4;
lattice_y = 1e4;

% Simulation parameters
% tau = 0.1; %base lag time
% time = 500; %simulation time in seconds
% time_pts = ceil(time/tau); %total time points
time_pts = 5000; %total time points (absolute time, camera frame-rate)
n = 25; %number of simulated particles.
random_start = 1; %0: all particles start at the center of the lattice, 1: particles are each assigned a random start location

%%% SET UP %%%
% Fetch a lattice
try
    disp('Attempting to load a pre-existing lattice.')
    load('lattice.mat','lattice')
    disp('Lattice loaded from file.')
catch
    disp('No pre-existing lattice is available. Generating a new lattice.')
    lattice = gen_lattice(save_lattice,heterogeneity,lattice_x,lattice_y);
    disp('Lattice generated successfully.')
end

if visualize_lattice == 1
    lattice_visualization = rescale(lattice,0,1);
    lattice_visualization = fliplr(rot90(lattice_visualization,-1));
    imshow(lattice_visualization)
end

% Failsafe in case a loaded lattice's settings don't match those of the current lattice size
lattice_x = size(lattice,2);
lattice_y = size(lattice,1);

% Convert lattice diffusivity values to (average) particle displacement values
displacements = round(100*sqrt(lattice)); %UNITS: 10^-2 um

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

for i = 1:n %iterate through each particle
    disp_message_particle = strcat(['Simulating particle #' num2str(i) '.']);
    disp(disp_message_particle)
    for j = 2:time_pts %for each particle, iterate through each time point
        try
            current_displacement_center = displacements(data_matrix(i,j-1,1),data_matrix(i,j-1,2));
        catch
            disp('WARNING. Particle struck the boundary and was rendered immobile.')
            data_matrix(i,j-1,1) = 0; %remove the displacement that crosses the boundary
            data_matrix(i,j-1,2) = 0; %remove the displacement that crosses the boundary
            break
        end
        
        current_displacement = current_displacement_center + round(get_dispmnt_variation(current_displacement_center,tau));
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
lattice_visualization = rescale(lattice,0,1);
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
    
    clearvars current_particle_x current_particle_y %clear variables for the next loop iteration34
end

%%% HISTOGRAM %%%
% Global displacement histogram plotting (all particles combined)
current_displacement_storage = reshape(current_displacement_storage,[1,numel(current_displacement_storage)]);
current_displacement_storage(current_displacement_storage==0) = [];

figure()
histogram(current_displacement_storage, 'Normalization', 'probability')
title('Step Size Distribution')
xlabel('Step Size')
ylabel('Frequency')

%%% MSD(DELTA-T) %%%
% MSD(delta-t) plot ** FIX FOR INCOMPLETE WALKS (BOUNDARY-STRUCK CASES)
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

%%% TIME-AVERAGED MSD(DELTA-TAU) %%%
% Time-averaged MSD(delta-tau) plot
% TBD
% multiples_delta_time = [1,5,10,50,100]; %time point intervals for displacement measurements
% displacement_jumps = zeros(n,time_pts,size(multiples_delta_time,2)); %storage for displacements at each (preset) time point interval