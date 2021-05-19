function lattice = gen_lattice(save_lattice,heterogeneity,lattice_size_x,lattice_size_y)
% GEN_LATTICE Generate a lattice environment for diffusion.
% save_lattice -- 0: don't save, 1: save (save the lattice to a .mat file)
% heterogeneity -- 0: perfectly homogeneous, 1: maximal heterogeneity
% lattice_size_x -- size of the lattice along the x-axis (width)
% lattice_size_y -- size of the lattice along the y-axis (height)

% Define the lattice
lattice_area = lattice_size_x * lattice_size_y;

% Initialize a void lattice
lattice = zeros(lattice_size_x,lattice_size_y);

% Get the number of diffusivity regions
num_regions = 8;
% num_regions = subregions(heterogeneity,lattice_area);

% Set diffusivity values for each region
% estimated via Wagner et al. Biomacromolecules article
min_diffusivity = 0.1; %um^2/s
max_diffusivity = 1.25; %um^2/s
% diffusivities = min_diffusivity + (max_diffusivity-min_diffusivity).*rand(1,num_regions);
diffusivities = linspace(min_diffusivity,max_diffusivity,num_regions); % ** MANUALLY FORCED HETEROGENEITY

% Generate the "root" points of each region. The regions will grow outward from these points.
region_start_pts = [randi([1,lattice_size_x],num_regions,1),randi([1,lattice_size_y],num_regions,1)];

% Grow the diffusivity regions until the entire lattice has been defined
for i = 1:lattice_size_x
    for j = 1:lattice_size_y
        dists_to_start_pts = pdist2(region_start_pts,[i,j]); %get the distance between each region start point and the current lattice index
        closest_dist = find(dists_to_start_pts==min(dists_to_start_pts),1); %find the closest region start point to the current lattice index
        lattice(i,j) = diffusivities(closest_dist); %set the diffusivity value for the given lattice index
    end
end

% Save the lattice
if save_lattice == 1
    save('lattice.mat','lattice')
end

end
