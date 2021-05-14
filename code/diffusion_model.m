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
n = 5; %number of simulated particles

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
displacements = ceil(rescale(sqrt(lattice),1,100));

% Matrix to store all relevant simulation data
data_matrix = zeros(n,time_pts,2); %for each particle at each time point, store the current x and y coordinates

% Set up particle start locations and initial diffusivities
data_matrix(:,1,1) = x_start;
data_matrix(:,1,2) = y_start;

for i = 1:n %iterate through each particle
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
            data_matrix(i,j:end,1) = data_matrix(i,j-1,1); % ** REMOVE THIS: KEEP FINAL VALUES AT ZERO AND REMOVE THEM PRIOR TO ANALYSIS (SHORTER TIMEFRAME)
            data_matrix(i,j:end,2) = data_matrix(i,j-1,2); % ** REMOVE THIS: KEEP FINAL VALUES AT ZERO AND REMOVE THEM PRIOR TO ANALYSIS (SHORTER TIMEFRAME)
            disp('WARNING. Particle struck the boundary and was rendered immobile.')
            break
        end
        
%         current_displacement = round(normrnd(current_displacement_center,10)); % **THE STANDARD DEVIATION IS ARBITRARY
        current_displacement = current_displacement_center + round(get_dispmnt_variation(current_displacement_center,tau));
        direction_select = randi(4); %randomly select a direction
            
%         temp_storage_actual(j-1) = current_displacement; %DEBUGGING
        
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
    plot(data_matrix(i,:,1),data_matrix(i,:,2)) %plot the entire trajectory
    plot(data_matrix(i,1,1),data_matrix(i,1,2),'>g','MarkerFaceColor','g','MarkerSize',10) %mark the start point
    plot(data_matrix(i,end,1),data_matrix(i,end,2),'sr','MarkerFaceColor','r','MarkerSize',10) %mark the end point
    legend('Particle Trajectory','Start','End')
    % xlim([0,lattice_y])
    % ylim([0,lattice_x])
    % xlabel('y')
    % ylabel('x')
end