% Test file

% test_lattice = zeros(3,3);
% 
% while ~all(all(test_lattice~=0))
%     for i = 1:3
%         for j = 1:3
%             test_lattice(i,j) = 1;
%         end
%     end
% end
% 
% test_lattice

%%%%%

% x = linspace(-100,100,10000);
% D = 100;
% tau = 0.1;
% 
% P = 1/(sqrt(4*pi*D*tau))*exp(-(x.^2/(4*D*tau)));
% 
% plot(x,P)
% title('PDF Given D=100, \tau=0.1')
% xlabel('x')
% ylabel('PDF')
% 
% cumtrapz(x,P)

%%%%%

% [x,pdf,cdf] = gen_PDF(100,0.1);
% plot(x,pdf)
% title('PDF')
% xlabel('Change to the step to be taken')
% ylabel('Probability')
% 
% figure()
% plot(x,cdf)
% 
% cdf = round(cdf,5);
% % figure()
% % plot(x,cdf)
% 
% displmnts = zeros(1000000,1);
% for i = 1:1000000
%     val = rand();
%     displmnt_idx = find(cdf==round(val,5),1);
%     displmnt = x(displmnt_idx);
%     try
%         displmnts(i)=displmnt;
%     catch
%         displmnts(i)=0;
%     end
% end
% 
% displmnts
% displmnts(displmnts==0) = [];
% 
% figure()
% histogram(displmnts,100)

%%%%%

% Time-averaged MSD(delta-tau) plot
% delta_taus = 1:(time_pts/10); %time point intervals for displacement measurements (given by the first 10% of time points)
% sqd_dispmnts_lag_time = zeros(n,time_pts,size(delta_taus,2)); %storage for displacements at each time point interval
% 
% % Compute the displacements for the given delta-tau values (multiples of delta-t)
% counter_msd_tau = 1;
% for dt = delta_taus
%     for i = 1:n
%         for j = 1:time_pts-dt
%             tamsd_plot_displacement_x = data_matrix(i,j+dt,1) - data_matrix(i,j,1);
%             tamsd_plot_displacement_y = data_matrix(i,j+dt,2) - data_matrix(i,j,2);
%             sqd_dispmnts_lag_time(i,j,counter_msd_tau) = tamsd_plot_displacement_x^2 + tamsd_plot_displacement_y^2;
%         end
%     end
%     counter_msd_tau = counter_msd_tau + 1;
% end
% 
% for i = 1:n
%     for j = 1:size(delta_taus,2)
%         if boundary_collision(i) == 1 %only modify the displacement data if the given  particle strikes the boundary
%             first_zero_idx = find(((sqd_dispmnts_lag_time(i,:,j)==0)+([diff(sqd_dispmnts_lag_time(i,:,j)) 0]==0))==2,1); %find the index of the two consecutive zeros in sqd_dispmnts_lag_time for the given multiple of delta-t
%             try
%                 sqd_dispmnts_lag_time(i,(first_zero_idx-delta_taus(j)):first_zero_idx-1,j) = 0; %remove the appropriate number of erroneous displacements
%             catch
%                 sqd_dispmnts_lag_time(i,:,j) = 0; %remove all displacements if the particle is immobilized at a time point lesser than the value of the time lag
%             end
%         end
%     end
% end

%%%%%

% a = [2,0,4,5,6,0,0,8,23,5,2,0,7,6,8,10,9,0,0,0,0,0,0,0,0,0,0,0];
% 
% N = 5; % Required number of consecutive numbers following a first one
% 
% x = diff(a)==0;
% 
% f = find([false,x]~=[x,false]);
% 
% g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first');
% 
% first_t = t(f(2*g-1)); % First t followed by >=N consecutive numbers

%%%%%

% test = [2,0,4,5,6,0,0,8,23,5,2,0,7,6,8,10,9,0,0,0,0,0,0,0,0,0,0,0];
% a = strfind(test,[0,0,0,0]);
% test(a(1):end) = [];

%%%%%

% D = 1000000;
% tau = 0.1;
% 
% % Define a sufficiently broad range of x values
% x = linspace(-1125,1125,1250000);
% 
% % Define the PDF given the diffusivity and time step values
% P = 1/(sqrt(4*pi*D*tau))*exp(-(x.^2/(4*D*tau)));
% 
% plot(x,P)
% title('PDF')
% xlabel('x')
% ylabel('PDF')

%%%%%

% function lattice = gen_lattice(save_lattice,heterogeneity,lattice_size_x,lattice_size_y)
% % GEN_LATTICE Generate a lattice environment for diffusion.
% % save_lattice -- 0: don't save, 1: save (save the lattice to a .mat file)
% % heterogeneity -- 0: perfectly homogeneous, 1: maximal heterogeneity
% % lattice_size_x -- size of the lattice along the x-axis (width)
% % lattice_size_y -- size of the lattice along the y-axis (height)
% 
% % Define the lattice
% lattice_area = lattice_size_x * lattice_size_y;
% 
% % Initialize a void lattice
% lattice = zeros(lattice_size_x,lattice_size_y);
% 
% % Get the number of diffusivity regions
% num_regions = 25;
% % num_regions = subregions(heterogeneity,lattice_area);
% 
% % Set diffusivity values for each region (estimated via Wagner et al. Biomacromolecules article)
% min_diffusivity = 0.1; %um^2/s
% max_diffusivity = 1.25; %um^2/s
% % diffusivities = min_diffusivity + (max_diffusivity-min_diffusivity).*rand(1,num_regions);
% diffusivities = linspace(min_diffusivity,max_diffusivity,num_regions); % ** MANUALLY FORCED HETEROGENEITY
% 
% % Generate the "root" points of each region. The regions will grow outward from these points.
% region_start_pts = [randi([1,lattice_size_x],num_regions,1),randi([1,lattice_size_y],num_regions,1)];
% 
% % Grow the diffusivity regions until the entire lattice has been defined
% for i = 1:lattice_size_x
%     for j = 1:lattice_size_y
%         dists_to_start_pts = pdist2(region_start_pts,[i,j]); %get the distance between each region start point and the current lattice index
%         closest_dist = find(dists_to_start_pts==min(dists_to_start_pts),1); %find the closest region start point to the current lattice index
%         lattice(i,j) = diffusivities(closest_dist); %set the diffusivity value for the given lattice index
%     end
%     
%     % Display lattice generation progress (every 10%)
%     if rem(i,(lattice_size_x/10))==0
%         lattice_disp_message = strcat(['Lattice generation: ' num2str(100*(i/lattice_size_x)) '% complete.']);
%         disp(lattice_disp_message)
%     end
% end
% 
% % Save the lattice
% if save_lattice == 1
%     save('lattice.mat','lattice')
%     save
% end
% 
% end

%%%%%

% D = 1000;
% tau = 0.1;
% x = linspace(-125,125,1250000);
% 
% cdf = gen_PDF(D,tau,x);

%%%%%
% num_regions = 25;
% % num_regions = subregions(heterogeneity,lattice_area);
% 
% % Set diffusivity values for each region (estimated via Wagner et al. Biomacromolecules article)
% min_diffusivity = 0.1; %um^2/s
% max_diffusivity = 1.25; %um^2/s
% 
% %Adjust the units of the diffusivities
% multiplier = 10000;
% min_diffusivity = multiplier*min_diffusivity; %10^-4 um^2/s
% max_diffusivity = multiplier*max_diffusivity; %10^-4 um^2/s
% 
% % diffusivities = min_diffusivity + (max_diffusivity-min_diffusivity).*rand(1,num_regions);
% diffusivities = round(linspace(min_diffusivity,max_diffusivity,num_regions)); % ** MANUALLY FORCED HETEROGENEITY
% 
% tau = 0.1;
%  
% x = linspace(-165,165,1650000);
% all_cdf = zeros(length(diffusivities),length(x));
% for i = 1:length(diffusivities)
%     all_cdf(i,:) = gen_PDF(diffusivities(i),tau,x);
% end
% 
% all_cdf(:,end)

%%%%%

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
    clearvars current_histogram_data_x current_histogram_data_y current_histogram_data_x_copy current_histogram_data_y_copy
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
        
        current_histogram_data_x = no_trailing_zeros(histogram_data_x(i,:,j));
        current_histogram_data_y = no_trailing_zeros(histogram_data_y(i,:,j));
        
        current_histogram_data_x_copy = no_trailing_zeros(histogram_data_x(i,:,j));
        current_histogram_data_y_copy = no_trailing_zeros(histogram_data_y(i,:,j));
        
%         idxs_zeros_x = find(~current_histogram_data_x);
%         idxs_zeros_y = find(~current_histogram_data_y);
        
        current_histogram_data_x((current_histogram_data_x_copy==0 & current_histogram_data_y_copy~=0)) = [];
        current_histogram_data_y((current_histogram_data_y_copy==0 & current_histogram_data_x_copy~=0)) = [];
        
        % Store the processed data for histogram plotting
        all_displacement_storage_x = [all_displacement_storage_x current_histogram_data_x];
        all_displacement_storage_y = [all_displacement_storage_y current_histogram_data_y];
        
    end
    
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