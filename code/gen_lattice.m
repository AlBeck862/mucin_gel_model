function [lattice,x,all_cdf,diffusivities] = gen_lattice(save_lattice,lattice_size_x,lattice_size_y,tau,single_diffusivity_toggle,single_diffusivity)
% GEN_LATTICE Generate a lattice environment for diffusion.
% save_lattice -- 0: don't save, 1: save (save the lattice to a .mat file)
% lattice_size_x -- size of the lattice along the x-axis (width)
% lattice_size_y -- size of the lattice along the y-axis (height)
% tau -- conversion factor: time points to seconds (unit of [seconds per time point])
% single_diffusivity -- 0: multiple subregions (heterogeneous), 1: uniform lattice (homogeneous)

% Define the lattice (homogeneous or heterogeneous)
if single_diffusivity_toggle == 1
    lattice = single_diffusivity*ones(lattice_size_x,lattice_size_y);
    diffusivities = single_diffusivity;
else
    % Initialize a void lattice
    lattice = zeros(lattice_size_x,lattice_size_y);

    % Get the number of diffusivity regions
    num_regions = 100;

    % Set diffusivity values for each region (estimated via Wagner et al. Biomacromolecules article)
%     min_diffusivity = 0.1; %um^2/s
%     max_diffusivity = 1.25; %um^2/s
    min_diffusivity = 0.05; %um^2/s
    max_diffusivity = 2.50; %um^2/s
    
    % Adjust the units of the diffusivities
    multiplier = 10000;
    min_diffusivity = multiplier*min_diffusivity; %10^-4 um^2/s
    max_diffusivity = multiplier*max_diffusivity; %10^-4 um^2/s

    % Compute the diffusivities of the lattice's subregions
    % diffusivities = min_diffusivity + (max_diffusivity-min_diffusivity).*rand(1,num_regions);
    diffusivities = round(linspace(min_diffusivity,max_diffusivity,num_regions)); % ** MANUALLY FORCED HETEROGENEITY
    
    % Generate the "root" points of each region. The regions will grow outward from these points.
    region_start_pts = [randi([1,lattice_size_x],num_regions,1),randi([1,lattice_size_y],num_regions,1)];

    % Grow the diffusivity regions until the entire lattice has been defined
    for i = 1:lattice_size_x
        for j = 1:lattice_size_y
            dists_to_start_pts = pdist2(region_start_pts,[i,j]); %get the distance between each region start point and the current lattice index
            closest_dist = find(dists_to_start_pts==min(dists_to_start_pts),1); %find the closest region start point to the current lattice index
            lattice(i,j) = diffusivities(closest_dist); %set the diffusivity value for the given lattice index
        end

        % Display lattice generation progress (every 10%)
        if rem(i,(lattice_size_x/10))==0
            lattice_disp_message = strcat(['Lattice generation: ' num2str(100*(i/lattice_size_x)) '% complete.']);
            disp(lattice_disp_message)
        end
    end
end

% Generate the CDF corresponding to each diffusivity value in the lattice
% x = linspace(-165,165,1650000);
x = linspace(-200,200,2000000);
all_cdf = zeros(length(diffusivities),length(x));
for i = 1:length(diffusivities)
    all_cdf(i,:) = gen_PDF(diffusivities(i),tau,x);
end

% Save the lattice
if save_lattice == 1
    save('lattice.mat','lattice')
    save('lattice_data.mat','x','all_cdf','diffusivities')
end

end