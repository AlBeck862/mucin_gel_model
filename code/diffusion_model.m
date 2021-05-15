% Diffusion model
% IMPORTANT NOTE: the x- and y-axis are reversed due to MATLAB row-column data storage convention

% Lattice parameters
visualize_lattice = 1; %0: no visualization, 1: visualization
save_lattice = 1; %0: don't save, 1: save (save the lattice to a .mat file if it is newly generated)
heterogeneity = 0.5; %0: perfectly homogeneous, 1: maximal heterogeneity **CURRENTLY NOT USED
lattice_x = 1e4;
lattice_y = 1e4;

% Simulation parameters
tau = 0.1; %base lag time
time = 500; %simulation time in seconds
time_pts = ceil(time/tau); %total time points
n = 3; %number of simulated particles.

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

% Particle start location **ADD A WAY TO RANDOMIZE THIS
x_start = round(lattice_x/2);
y_start = round(lattice_y/2);

% Rescale lattice for motion compatibility
displacements = round(100*sqrt(lattice)); %UNITS: 10^-2 um

% displacements(8000:8025,1:10)
% min(min(displacements))
% max(max(displacements))
% dsbfiuas
% Matrix to store all relevant simulation data
data_matrix = zeros(n,time_pts,2); %for each particle at each time point, store the current x and y coordinates

% Set up particle start locations and initial diffusivities
data_matrix(:,1,1) = x_start;
data_matrix(:,1,2) = y_start;

% Storage for the global histogram of displacements
current_displacement_storage = zeros(n,time_pts-1);

for i = 1:n %iterate through each particle
    disp_message_particle = strcat(['Simulating particle #' num2str(i) '.']);
    disp(disp_message_particle)
%     temp_storage_actual = zeros(1,time_pts-1); %DEBUGGING
%     temp_storage_base = zeros(1,time_pts-1); %DEBUGGING
%     temp_storage_lattice = zeros(1,time_pts-1); %DEBUGGING
%     temp_storage_x = zeros(1,time_pts-1); %DEBUGGING
%     temp_storage_y = zeros(1,time_pts-1); %DEBUGGING
    for j = 2:time_pts %for each particle, iterate through each time point
        try
            current_displacement_center = displacements(data_matrix(i,j-1,1),data_matrix(i,j-1,2));
%             temp_storage_base(j-1) = current_displacement_center; %DEBUGGING
%             temp_storage_lattice(j-1) = lattice(data_matrix(i,j-1,1),data_matrix(i,j-1,2)); %DEBUGGING
%             temp_storage_x(j-1) = data_matrix(i,j-1,1); %DEBUGGING
%             temp_storage_y(j-1) = data_matrix(i,j-1,2); %DEBUGGING
        catch
%             data_matrix(i,j:end,1) = data_matrix(i,j-1,1); % ** REMOVE THIS: KEEP FINAL VALUES AT ZERO AND REMOVE THEM PRIOR TO ANALYSIS (SHORTER TIMEFRAME)
%             data_matrix(i,j:end,2) = data_matrix(i,j-1,2); % ** REMOVE THIS: KEEP FINAL VALUES AT ZERO AND REMOVE THEM PRIOR TO ANALYSIS (SHORTER TIMEFRAME)
            disp('WARNING. Particle struck the boundary and was rendered immobile.')
            break
        end
        
        current_displacement = current_displacement_center + round(get_dispmnt_variation(current_displacement_center,tau));
        direction_select = randi(4); %randomly select a direction
            
        current_displacement_storage(i,j-1) = current_displacement;
        
        if direction_select == 1 %+x (RIGHT)
            data_matrix(i,j,1) = data_matrix(i,j-1,1) + abs(current_displacement);
            data_matrix(i,j,2) = data_matrix(i,j-1,2);
        elseif direction_select == 2 %-x (LEFT)
            data_matrix(i,j,1) = data_matrix(i,j-1,1) - abs(current_displacement);
            data_matrix(i,j,2) = data_matrix(i,j-1,2);
        elseif direction_select == 3 %+y (DOWN)
            data_matrix(i,j,1) = data_matrix(i,j-1,1);
            data_matrix(i,j,2) = data_matrix(i,j-1,2) + abs(current_displacement);
        elseif direction_select == 4 %-y (UP)
            data_matrix(i,j,1) = data_matrix(i,j-1,1);
            data_matrix(i,j,2) = data_matrix(i,j-1,2) - abs(current_displacement);
        end
    end
end

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
    % xlim([0,lattice_y])
    % ylim([0,lattice_x])
    % xlabel('y')
    % ylabel('x')
    clearvars current_particle_x current_particle_y %clear variables for the next loop iteration34
end

% Global displacement histogram plotting (all particles combined)
current_displacement_storage = reshape(current_displacement_storage,[1,numel(current_displacement_storage)]);
current_displacement_storage(current_displacement_storage==0) = [];

figure()
histogram(current_displacement_storage)
title('Step Size Distribution')
xlabel('Step Size')
ylabel('Frequency')

%%% WORK IN PROGRESS %%% 
% % Processing loop (all data analysis goes inside this loop)
% particle_data_x = 0; %dummy initialization to allow for the clearvars statement to work on the first loop iteration
% particle_data_y = 0; %dummy initialization to allow for the clearvars statement to work on the first loop iteration
% for i = 1:n
%     clearvars particle_data_x particle_data_y %clear the particle_data variable for each new particle
% 
%     particle_data_x = data_matrix(i,:,1); %storage the x coordinates
%     particle_data_y = data_matrix(i,:,2); %storage the y coordinates
%     
%     particle_data_x(particle_data_x==0) = []; %remove extra zeros in case the particle was immobilized by the boundary
%     particle_data_y(particle_data_y==0) = []; %remove extra zeros in case the particle was immobilized by the boundary
%     
%     %DO WE ONLY CARE ABOUT DISPLACEMENT IN X?? OR IN BOTH X AND Y?
%     
%     
% end