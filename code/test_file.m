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

% %%% HISTOGRAMS %%%
% % Histogram for all displacements (1*(delta-t))
% clearvars all_displacement_storage
% all_displacement_storage = [];
% for i = 1:n
%     % Store the processed data for histogram plotting
%     all_displacement_storage = [all_displacement_storage no_trailing_zeros(current_displacement_storage(i,:))];
% end
% 
% figure()
% histogram(all_displacement_storage, 'Normalization', 'pdf')
% title('Step Size Distribution')
% xlabel('\Deltax, \Deltay [10^{-2}\mum]')
% ylabel('P(\Deltax, \Deltay, \Delta\tau=1 time point)')
% 
% % Histograms for displacements using greater multiples of delta-t (iteratively: multiples_delta_time*(delta-t))
% histogram_data_x = zeros(n,time_pts,size(multiples_delta_time,2));
% histogram_data_y = zeros(n,time_pts,size(multiples_delta_time,2));
% counter_hist = 1;
% for dt = multiples_delta_time
%     for i = 1:n
%         for j = 1:time_pts-dt
%             histogram_data_x(i,j,counter_hist) = data_matrix(i,j+dt,1) - data_matrix(i,j,1);
%             histogram_data_y(i,j,counter_hist) = data_matrix(i,j+dt,2) - data_matrix(i,j,2);
%         end
%     end
%     counter_hist = counter_hist + 1;
% end
% 
% for j = 1:size(multiples_delta_time,2)
%     clearvars all_displacement_storage_x all_displacement_storage_y all_displacement_storage_both
%     clearvars current_histogram_data_x current_histogram_data_y current_histogram_data_x_copy current_histogram_data_y_copy
%     all_displacement_storage_x = [];
%     all_displacement_storage_y = [];
%     for i = 1:n
%         % Remove erroneous displacements that result from a particle striking the boundary
%         if boundary_collision(i) == 1 %only modify the displacement data if the given particle strikes the boundary
%             start_of_trailing_zeros = find(histogram_data_x(i,:,j),1,'last') + 1;
%             try
%                 % Remove the appropriate number of erroneous displacements
%                 histogram_data_x(i,(start_of_trailing_zeros-multiples_delta_time(j)):start_of_trailing_zeros-1,j) = 0; 
%                 histogram_data_y(i,(start_of_trailing_zeros-multiples_delta_time(j)):start_of_trailing_zeros-1,j) = 0;
%             catch
%                 % If the particle stopped so early that all displacements are erroneous, set the enter set of displacements to zero
%                 histogram_data_x(i,:,j) = 0;
%                 histogram_data_y(i,:,j) = 0;
%             end
%         end
%         
%         current_histogram_data_x = no_trailing_zeros(histogram_data_x(i,:,j));
%         current_histogram_data_y = no_trailing_zeros(histogram_data_y(i,:,j));
%         
%         current_histogram_data_x_copy = no_trailing_zeros(histogram_data_x(i,:,j));
%         current_histogram_data_y_copy = no_trailing_zeros(histogram_data_y(i,:,j));
%         
% %         idxs_zeros_x = find(~current_histogram_data_x);
% %         idxs_zeros_y = find(~current_histogram_data_y);
%         
%         current_histogram_data_x((current_histogram_data_x_copy==0 & current_histogram_data_y_copy~=0)) = [];
%         current_histogram_data_y((current_histogram_data_y_copy==0 & current_histogram_data_x_copy~=0)) = [];
%         
%         % Store the processed data for histogram plotting
%         all_displacement_storage_x = [all_displacement_storage_x current_histogram_data_x];
%         all_displacement_storage_y = [all_displacement_storage_y current_histogram_data_y];
%         
%     end
%     
%     all_displacement_storage_both = [all_displacement_storage_x all_displacement_storage_y];
%     
% %     figure()
% %     histogram(all_displacement_storage_x, 'Normalization', 'pdf')
% %     title('Step Size Distribution')
% %     xlabel('\Deltax [10^{-2}\mum]')
% %     hist_y_label_str = strcat(['P(\Deltax, \Delta\tau=' num2str(multiples_delta_time(j)) ' time points)']);
% %     ylabel(hist_y_label_str)
% %     
% %     figure()
% %     histogram(all_displacement_storage_y, 'Normalization', 'pdf')
% %     title('Step Size Distribution')
% %     xlabel('\Deltay [10^{-2}\mum]')
% %     hist_y_label_str = strcat(['P(\Deltay, \Delta\tau=' num2str(multiples_delta_time(j)) ' time points)']);
% %     ylabel(hist_y_label_str)
% 
% 	figure()
%     histogram(all_displacement_storage_both, 'Normalization', 'pdf')
%     title('Step Size Distribution')
%     xlabel('\Deltax, \Deltay [10^{-2}\mum]')
%     hist_y_label_str = strcat(['P(\Deltax, \Deltay, \Delta\tau=' num2str(multiples_delta_time(j)) ' time points)']);
%     ylabel(hist_y_label_str)
% end

%%%%%

% save_lattice = 0;
% heterogeneity = 1;
% lattice_size_x = 1000;
% lattice_size_y = 1000;
% tau = 0.1;
% tic
% [lattice,x,all_cdf,diffusivities] = gen_lattice(save_lattice,heterogeneity,lattice_size_x,lattice_size_y,tau);
% toc
% 
% % Time without parallel computing: 27.44 seconds / 26.90 seconds
% % Time with parallel computing: 49.93 seconds / 38.92 seconds

%%%%%

% time_idx = [1,2,3,4,5,10,50,100,125,150,175,200];
% data_10000 = [0.996,0.8192,0.7588,0.7343,0.723,0.713,0.712,0.71,0.714,0.718,0.723,0.726];
% data_2500 = [0.992,0.822,0.762,0.736,0.723,0.711,0.718,0.707,0.702,0.699,0.696,0.692];
% data_2500_again = [0.996,0.823,0.764,0.739,0.726,0.714,0.716,0.718,0.724,0.727,0.728,0.729];
% 
% raw_theory_10000 = [44.72,63.25,77.46,89.44,100.00,141.42,316.23,447.21,500,547.72,591.61,632.46];
% raw_theory_2500 = [22.36,31.62,38.73,44.72,50.00,70.71,158.11,223.61,250.00,273.86,295.80,316.23];
% raw_10000 = [44.54,51.82,58.78,65.68,72.26,100.90,225.14,317.58,356.79,393.52,427.60,459.20];
% raw_2500 = [22.19,25.98,29.53,32.90,36.17,50.27,113.49,158.01,175.46,191.38,205.89,218.82];
% raw_2500_again = [22.26,26.02,29.59,33.05,36.31,50.50,113.23,160.59,180.91,199.02,215.38,230.49];
% 
% raw_2500_again_again = [15.816,22.3997,27.4605,31.7248,35.4858,50.3993,115.0598,163.575,182.434,199.6032,215.3026,229.548];

% Unchanged
% figure()
% plot(time_idx,data_10000,'-o')
% hold on
% plot(time_idx,data_2500,'-o')
% plot(time_idx,data_2500_again,'-o')
% legend('D=10000','D=2500','D=2500 again')
% xlabel('# time points')
% ylabel('Histogram std. dev. / Theoretical std. dev.')
% 
% figure()
% plot(time_idx,raw_theory_10000,'--')
% hold on
% plot(time_idx,raw_theory_2500,'--')
% plot(time_idx,raw_10000,'-o')
% plot(time_idx,raw_2500,'-o')
% plot(time_idx,raw_2500_again,'-o')
% plot(time_idx,raw_2500_again_again,'-o')
% legend('Theory D=10000','Theory D=2500','D=10000','D=2500','D=2500 again','D=2500 again again')
% xlabel('# time points')
% ylabel('Std. Deviation')
% 
% figure()
% plot(raw_theory_10000,raw_10000,'-o')
% hold on
% plot(raw_theory_2500,raw_2500,'-o')
% plot(raw_theory_2500,raw_2500_again,'-o')
% plot(raw_theory_2500,raw_2500_again_again,'-o')
% plot(1:500,'--')
% legend('D=10000','D=2500','D=2500 again','D=2500 again again','Perfect Agreement')
% xlabel('Theoretical Std. Dev.')
% ylabel('Experimental Std. Dev.')


% New test
% raw_2500_first_method = [22.2205,26.0039,29.5289,32.9218,36.1586,50.1428,111.8225,154.4281,171.5639,187.2045,201.7564,214.9222];
% raw_2500_second_method = [15.8517,22.4354,27.4623,31.6767,35.3775,49.9109,111.6039,154.251,171.3927,186.9898,201.5469,214.7159];
% raw_2500_third_method = [15.8481,22.4294,27.4527,31.6631,35.3591,49.861,111.0449,152.7034,169.2417,184.1718,198.0009,210.3952];
% 
% figure()
% plot(raw_theory_2500,raw_2500_first_method,'-o')
% hold on
% plot(raw_theory_2500,raw_2500_second_method,'-o')
% plot(raw_theory_2500,raw_2500_third_method,'-o')
% plot(1:500,'--')
% legend('First Method','Second Method','Third Method','Perfect Agreement')
% xlabel('Theoretical Std. Dev.')
% ylabel('Experimental Std. Dev.')
% 
% fit_line = polyfit(raw_theory_2500,raw_2500_third_method,1)
% quotients = raw_2500_third_method./raw_theory_2500
% avg_quotient = mean(quotients)

% Squared
% figure()
% plot(time_idx,raw_theory_10000.^2,'--')
% hold on
% plot(time_idx,raw_theory_2500.^2,'--')
% plot(time_idx,raw_10000.^2,'-o')
% plot(time_idx,raw_2500.^2,'-o')
% plot(time_idx,raw_2500_again.^2,'-o')
% legend('Theory D=10000','Theory D=2500','D=10000','D=2500','D=2500 again')
% xlabel('# time points')
% ylabel('Std. Deviation')
% 
% figure()
% plot(raw_theory_10000.^2,raw_10000.^2,'-o')
% hold on
% plot(raw_theory_2500.^2,raw_2500.^2,'-o')
% plot(raw_theory_2500.^2,raw_2500_again.^2,'-o')
% plot(1:500000,'--')
% legend('D=10000','D=2500','D=2500 again','Perfect Agreement')
% xlabel('Theoretical Std. Dev.')
% ylabel('Experimental Std. Dev.')

% Fitting
% sqrt_fit_raw_theory_10000 = polyfit(sqrt(time_idx),raw_theory_10000,1)
% sqrt_fit_raw_10000 = polyfit(sqrt(time_idx),raw_10000,1)
% 
% fit_raw_10000 = polyfit(raw_theory_10000.^2,raw_10000.^2,1)

%%%%%

% try
%     delete(pool)
% end
% 
% num = 3; %number of total iterations
% workers = 3; %number of cores
% pool = parpool(workers);
% result = zeros(num,2);
% 
% for round = 1:ceil(num/workers)
%     if round*workers<=num
%         val = workers;
%         spmd (val)
%             a = zeros(1,2);
%             disp(labindex)
%             a(1) = round*labindex;
%             a(2) = round*labindex;
%         end
%     else
%         val = rem(num,workers);
%         spmd (val)
%             disp(labindex)
%             a = round*labindex;
%         end
%     end
%     
% %     [a{:}]
%     idx = find(~result,1,'first');
%     for lol = 1:val
%         result(idx+lol-1,:) = [a{lol}];
%     end
% %     clear a
%     
% end

%%%%%

num_regions = 100;
tau = 0.1;

% Set diffusivity values for each region (estimated via Wagner et al. Biomacromolecules article)
%     min_diffusivity = 0.1; %um^2/s
%     max_diffusivity = 1.25; %um^2/s
min_diffusivity = 0.05; %um^2/s
max_diffusivity = 3.00; %um^2/s

% Adjust the units of the diffusivities
multiplier = 10000;
min_diffusivity = multiplier*min_diffusivity; %10^-4 um^2/s
max_diffusivity = multiplier*max_diffusivity; %10^-4 um^2/s

% Compute the diffusivities of the lattice's subregions
% diffusivities = min_diffusivity + (max_diffusivity-min_diffusivity).*rand(1,num_regions);
diffusivities = round(linspace(min_diffusivity,max_diffusivity,num_regions)); % ** MANUALLY FORCED HETEROGENEITY

x = linspace(-200,200,2000000);
all_cdf = zeros(length(diffusivities),length(x));
for i = 1:length(diffusivities)
    all_cdf(i,:) = gen_PDF(diffusivities(i),tau,x);
end
