tic %begin benchmarking (entire script)

%%% PARAMETERS %%%
% Lattice generation parameters
% Note: the origin (i.e. (0,0)) of all images is in the top left corner due to MATLAB's data storage convention
lattice_x = 7.5e4;                            %size of the lattice along the horizontal axis (number of lattice columns)
lattice_y = 7.5e4;                            %size of the lattice along the vertical axis (number of lattice rows)
save_lattice = true;                        %false: don't save, true: save (save the lattice to a .mat file if it is newly generated)
single_diffusivity_toggle = false;          %false: multiple subregions, true: uniform lattice (note: setting "import_lattice" to "true" bypasses this setting)
single_diffusivity = 1000;                 %value of the diffusivity when constructing a single-diffusivity lattice (units: 10^-4 um^2/s)
import_lattice = true;                      %false: automatic lattice generation, true: import a manually-designed lattice
invert_grayscale = true;                   %false: do not invert the grayscale image, true: invert the grayscale image
round_imported_lattice = true;              %false: do not modify the grayscale image's pixel values, true: round the grayscale image's pixel values to the nearest multiple of "round_imported_lattice_multiple"
round_imported_lattice_multiple = 5;        %round the grayscale image's pixel values to the nearest multiple of this value
multiplier = 10000;                           %lattice scaling: multiply the provided units (e.g. micrometer) by the inverse of this value

% Simulation parameters
time_pts = 5000;                            %total time points (absolute time, camera frame-rate)
n = 25;                                      %number of simulated particles.
start_type = 'random';                      %'center': all particles start at the center of the lattice, 'random': particles are each assigned a random start location, 'other_fixed': all particles start at a hard-coded location
visualize_lattice = true;                   %false: no visualization, true: visualization
conversion_factor = 0.1;                    %conversion factor, units of seconds per time point

% Plotting parameters
multiples_delta_time = [1,10,100];          %additional time point intervals for displacement histogram generation (each value corresponds to a set of two histograms (x-direction and y-direction)
msd_dtau_log = true;                        %false: MSD(delta-tau) plots will be linearly scaled, true: MSD(delta-tau) plots will be logarithmically scaled (this does not affect the line-of-best-fit parameters)
moving_avg_kernel = round(time_pts/100);    %size of the window over which to average the displacement data when plotting step-size over time

% Other parameters
save_data = false;                          %false: no data variables are saved, true: certain data variables are saved (data_matrix, boundary_collision, all_displacement_storage_x, all_displacement_storage_y)

%%% SET-UP %%%
% Fetch a lattice
try
    disp('Attempting to load a pre-existing lattice and associated data.')
    load('lattice.mat','lattice')
    load('lattice_data.mat','x','all_cdf','diffusivities')
    disp('Lattice and associated data loaded from file.')
catch
    tic %begin benchmarking (lattice generation)
    [lattice,x,all_cdf,diffusivities] = gen_lattice(save_lattice,lattice_x,lattice_y,conversion_factor,single_diffusivity_toggle,single_diffusivity,import_lattice, invert_grayscale, round_imported_lattice, round_imported_lattice_multiple, multiplier);
    toc %end benchmarking (lattice generation)
end

%%% LATTICE VISUALIZATION %%%
% Display and save the lattice visualization if the image is small enough
if lattice_x*lattice_y <= 1e8
    lattice_visualization(visualize_lattice,lattice)
end

%%% WALK SIMULATION %%%
[data_matrix,boundary_collision] = walk_simulation(n,time_pts,start_type,lattice,save_data,x,all_cdf,diffusivities);

%%% WALK VISUALIZATION %%%
% Display and save the walk visualizations if the image is small enough
if lattice_x*lattice_y <= 1e8
    walk_visualization(n,lattice,data_matrix)
end

%%% HISTOGRAMS %%%
histogram_plotting(n,time_pts,multiples_delta_time,data_matrix,boundary_collision,save_data,conversion_factor,multiplier)

%%% MSD(DELTA-T) %%%
msd_dt_plotting(n,time_pts,data_matrix,conversion_factor,multiplier)

%%% TIME-AVERAGED MSD(DELTA-TAU) %%%
msd_dtau_plotting(n,time_pts,data_matrix,boundary_collision,conversion_factor,msd_dtau_log,multiplier)

%%% STEP-SIZE OVER TIME & DIFFUSIVITY FREQUENCIES %%%
stepsize_dt_plotting(n,time_pts,data_matrix,lattice,moving_avg_kernel,conversion_factor,multiplier,diffusivities)

toc %end benchmarking (entire script)