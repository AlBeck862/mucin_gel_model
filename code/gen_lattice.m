function [lattice,diffusivities] = gen_lattice(save_lattice,lattice_size_x,lattice_size_y,tau,single_diffusivity_toggle,single_diffusivity,import_lattice, invert_grayscale, round_imported_lattice, round_imported_lattice_multiple, multiplier)
% GEN_LATTICE Generate a lattice environment for diffusion.

%%% DIFFUSIVITY PARAMETERS %%%
% % min_diffusivity = 0.1;  %units: um^2/s
% % max_diffusivity = 1.25; %units: um^2/s

% Particle radius: r=60nm
min_diffusivity = 0.041;  %units: um^2/s
max_diffusivity = 4.090; %units: um^2/s

% Particle radius: r=500nm
% min_diffusivity = 0.005;  %units: um^2/s
% max_diffusivity = 0.490; %units: um^2/s

% Adjust the units of the diffusivities
min_diffusivity = multiplier*min_diffusivity; %units: 10^-4 um^2/s
max_diffusivity = multiplier*max_diffusivity; %units: 10^-4 um^2/s

%%% LATTICE IMPORT OR GENERATION %%%
if import_lattice == false %automatic lattice generation
	disp('No pre-existing lattice is available. Generating a new lattice.')
    
    % Define the lattice (homogeneous or heterogeneous)
    if single_diffusivity_toggle == true
        lattice = single_diffusivity*ones(lattice_size_x,lattice_size_y);
        diffusivities = single_diffusivity;
    else
        % Initialize a void lattice
        lattice = zeros(lattice_size_x,lattice_size_y);

        % Get the number of diffusivity regions
        num_regions = 25;

        % Compute the diffusivities of the lattice's subregions
        diffusivities = round(linspace(min_diffusivity,max_diffusivity,num_regions));

        % Generate the "root" points of each region. The regions will grow outward from these points.
        region_start_pts = [randi([1,lattice_size_x],num_regions,1),randi([1,lattice_size_y],num_regions,1)];

        % Grow the diffusivity regions until the entire lattice has been defined
        for i = 1:lattice_size_x
            for j = 1:lattice_size_y
                dists_to_start_pts = pdist2(region_start_pts,[i,j]);                %get the distance between each region start point and the current lattice index
                closest_dist = find(dists_to_start_pts==min(dists_to_start_pts),1); %find the closest region start point to the current lattice index
                lattice(i,j) = diffusivities(closest_dist);                         %set the diffusivity value for the given lattice index
            end

            % Display lattice generation progress (every 10%)
            if rem(i,(lattice_size_x/10))==0
                lattice_disp_message = strcat(['Lattice generation: ' num2str(100*(i/lattice_size_x)) '% complete.']);
                disp(lattice_disp_message)
            end
        end
    end
	disp('Lattice generated successfully.')

else %manually-designed lattice import
    disp('No pre-existing lattice is available. Importing the lattice from the provided image.')
    
    % Import the source image and convert it to grayscale
    lattice_img = imread('lattice.png');        %load the source image
	gray_image = im2gray(lattice_img);          %convert the image to grayscale
    gray_image = fliplr(rot90(gray_image,-1));  %ensure the lattice will retain the source image's orientation

    % Resize the image
    original_sz = size(gray_image);                                 %get the size of the original image
    expansion_factor_x = lattice_size_x/original_sz(1);             %factor by which the image will be expanded along the x-axis
    expansion_factor_y = lattice_size_y/original_sz(2);             %factor by which the image will be expanded along the y-axis
    xg = 1:original_sz(1);                                          %define the x-direction size of the original image for the interpolation function
    yg = 1:original_sz(2);                                          %define the y-direction size of the original image for the interpolation function
    F = griddedInterpolant({xg,yg},double(gray_image),'nearest');   %define the interpolation function
    xq = (0:1/expansion_factor_x:original_sz(1))';                  %define the x-direction size of the resized image for the interpolation function
    yq = (0:1/expansion_factor_y:original_sz(2))';                  %define the x-direction size of the resized image for the interpolation function
    vq = uint8(F({xq,yq}));                                         %resize the image
    lattice = vq(1:end-1,1:end-1);                                  %trim the resized image to the correct dimensions
    final_sz = size(lattice);                                       %get the size of the processed image
    
    whos lattice
    
    % Invert the grayscale image if requested
    if invert_grayscale == true
        lattice = imcomplement(lattice);
    end
    
    % Round the grayscale image's pixel values if requested
    if round_imported_lattice == true
%         lattice = single(lattice);
        lattice = round(lattice/round_imported_lattice_multiple)*round_imported_lattice_multiple; %round to the nearest multiple of "round_imported_lattice_multiple"
    end
    
    lattice_disp_message = strcat(['Successfully imported the lattice and converted it from a ' num2str(original_sz(1)) '-by-' num2str(original_sz(2)) ' matrix to a ' num2str(final_sz(1)) '-by-' num2str(final_sz(2)) ' matrix.']);
    disp(lattice_disp_message)
    
    % Convert grayscale values to diffusivities
    lattice = uint16(round(rescale(lattice,min_diffusivity,max_diffusivity)));
    
    whos lattice
    
    diffusivities = double(unique(lattice));
end

% Generate the CDF corresponding to each diffusivity value in the lattice
% x = linspace(-300,300,3000000);
% all_cdf = zeros(length(diffusivities),length(x));
% for i = 1:length(diffusivities)
%     all_cdf(i,:) = gen_PDF(diffusivities(i),tau,x);
% end

% Save the lattice
if save_lattice == true
    save('lattice.mat','lattice')
    save('diffusivities.mat','diffusivities')
end

end