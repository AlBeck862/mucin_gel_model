tic %begin benchmarking (entire script)

% Diffusion model
% IMPORTANT NOTE: the x-axis and y-axis are reversed due to MATLAB's row-column data storage convention

global x all_cdf diffusivities

%%% PARAMETERS %%%
% Lattice generation parameters
save_lattice = 1; %0: don't save, 1: save (save the lattice to a .mat file if it is newly generated) --> used only if a lattice is generated
lattice_x = 1e4; %size of the lattice along the horizontal axis (number of lattice columns)
lattice_y = 1e4; %size of the lattice along the vertical axis (number of lattice rows)
single_diffusivity_toggle = 1; %0: multiple subregions, 1: uniform lattice
single_diffusivity = 1000; %value of the diffusivity when constructing a single-diffusivity lattice

% Simulation parameters
time_pts = 20000; %total time points (absolute time, camera frame-rate)
n = 5; %number of simulated particles.
random_start = 0; %0: all particles start at the center of the lattice, 1: particles are each assigned a random start location, -1: all particles start at a hard-coded location
visualize_lattice = 1; %0: no visualization, 1: visualization
conversion_factor = 0.1; %conversion factor, units of seconds per time point

% Plotting parameters
% multiples_delta_time = [1,5,10,50,100,150,200]; %additional time point intervals for displacement histogram generation (each value corresponds to a set of two histograms (x-direction and y-direction)
multiples_delta_time = [1,10,100]; %additional time point intervals for displacement histogram generation (each value corresponds to a set of two histograms (x-direction and y-direction)
msd_dtau_log = 1; %0: MSD(delta-tau) plots will be linearly scaled, 1: MSD(delta-tau) plots will be logarithmically scaled (this does not affect the line-of-best-fit parameters)
moving_avg_kernel = round(time_pts/100); %size of the window over which to average the displacement data when plotting step-size over time

% Other parameters
save_data = 0; %0: no data variables are saved, 1: certain data variables are saved (data_matrix, boundary_collision, all_displacement_storage_x, all_displacement_storage_y)

%%% SET-UP %%%
% Fetch a lattice
try
    disp('Attempting to load a pre-existing lattice and associated data.')
    load('lattice.mat','lattice')
    load('lattice_data.mat','x','all_cdf','diffusivities')
    disp('Lattice and associated data loaded from file.')
catch
    disp('No pre-existing lattice is available. Generating a new lattice.')
    tic %begin benchmarking (lattice generation)
    [lattice,x,all_cdf,diffusivities] = gen_lattice(save_lattice,lattice_x,lattice_y,conversion_factor,single_diffusivity_toggle,single_diffusivity);
    toc %end benchmarking (lattice generation)
    disp('Lattice generated successfully.')
end

%%% LATTICE VISUALIZATION %%%
% lattice_visualization(visualize_lattice,lattice)

%%% WALK SIMULATION %%%
% Simulate the random walk
[data_matrix,boundary_collision] = walk_simulation(n,time_pts,random_start,lattice,save_data);

%%% WALK VISUALIZATION %%%
% walk_visualization(n,lattice,data_matrix)

%%% HISTOGRAMS %%%
histogram_plotting(n,time_pts,multiples_delta_time,data_matrix,boundary_collision,save_data)

%%% MSD(DELTA-T) %%%
msd_dt_plotting(n,time_pts,data_matrix)

%%% TIME-AVERAGED MSD(DELTA-TAU) %%%
msd_dtau_plotting(n,time_pts,data_matrix,boundary_collision,conversion_factor,msd_dtau_log)

%%% STEP-SIZE OVER TIME & DIFFUSIVITY FREQUENCIES %%%
% stepsize_dt_plotting(n,data_matrix,lattice,moving_avg_kernel)

toc %end benchmarking (entire script)