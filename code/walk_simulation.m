function [data_matrix,boundary_collision] = walk_simulation(n,time_pts,start_type,lattice,save_data,x,all_cdf,diffusivities)
% WALK_SIMULATION Simulate the diffusion of n particles over a maximum of
% time_pts time points.
% start_type -- determines the initial location of the particle (see the main script for details)
% lattice -- stores the diffusivity values used to define a particle's diffusion given its location at a given time point

tic %begin benchmarking (particle simulation)

% Failsafe in case a loaded lattice's settings don't match those of the current lattice size
lattice_x = size(lattice,1);
lattice_y = size(lattice,2);

% Matrix to store all relevant simulation data
data_matrix = zeros(n,time_pts,2); %for each particle at each time point, store the current x and y coordinates

% Memory of boundary collisions
boundary_collision = zeros(1,n);

for i = 1:n %iterate through each particle
    disp_message_particle = strcat(['Simulating particle #' num2str(i) '.']);
    disp(disp_message_particle)
    
    data_matrix(i,1,:) = start_location(start_type,lattice_x,lattice_y); %set the initial position of the particle
    
    for j = 2:time_pts %for each particle, iterate through each time point
        try
            current_diffusivity = lattice(data_matrix(i,j-1,1),data_matrix(i,j-1,2));
        catch
            disp('WARNING. The particle struck the boundary and was rendered immobile.')
            data_matrix(i,j-1,1) = 0;   %remove the displacement that crosses the boundary
            data_matrix(i,j-1,2) = 0;   %remove the displacement that crosses the boundary
            boundary_collision(i) = 1;  %remember that this particle struck the boundary
            break
        end
        
%         tic
        current_displacement_x = round(get_dispmnt_variation(current_diffusivity,x,all_cdf,diffusivities)); %randomly select a distance in the x-direction
        current_displacement_y = round(get_dispmnt_variation(current_diffusivity,x,all_cdf,diffusivities)); %randomly select a distance in the y-direction
%         toc
        
        % Update the particle's position
        data_matrix(i,j,1) = data_matrix(i,j-1,1) + current_displacement_x;
        data_matrix(i,j,2) = data_matrix(i,j-1,2) + current_displacement_y;
    end
end

if save_data == true
    save('walk_data.mat','data_matrix','boundary_collision')
end

toc %end benchmarking (particle simulation)

end