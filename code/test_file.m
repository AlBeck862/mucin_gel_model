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

% num_regions = 100;
% tau = 0.1;
% 
% % Set diffusivity values for each region (estimated via Wagner et al. Biomacromolecules article)
% %     min_diffusivity = 0.1; %um^2/s
% %     max_diffusivity = 1.25; %um^2/s
% min_diffusivity = 0.05; %um^2/s
% max_diffusivity = 3.00; %um^2/s
% 
% % Adjust the units of the diffusivities
% multiplier = 10000;
% min_diffusivity = multiplier*min_diffusivity; %10^-4 um^2/s
% max_diffusivity = multiplier*max_diffusivity; %10^-4 um^2/s
% 
% % Compute the diffusivities of the lattice's subregions
% % diffusivities = min_diffusivity + (max_diffusivity-min_diffusivity).*rand(1,num_regions);
% diffusivities = round(linspace(min_diffusivity,max_diffusivity,num_regions)); % ** MANUALLY FORCED HETEROGENEITY
% 
% x = linspace(-200,200,2000000);
% all_cdf = zeros(length(diffusivities),length(x));
% for i = 1:length(diffusivities)
%     all_cdf(i,:) = gen_PDF(diffusivities(i),tau,x);
% end

%%%%%

% fh = figure;
% imageh = imshow(false(1000));
% 
% % Create buttons in the figure.
% uicontrol('Parent',fh,'Style','pushbutton','String','paint','Callback',{@paintButtonCallback, imageh});
% bh = uicontrol('Parent',fh,'Style','pushbutton','String','line','Callback',{@lineButtonCallback, imageh});
% bh.Position(2) = 50;
% bh2 = uicontrol('Parent',fh,'Style','pushbutton','String','line2','Callback',{@line2ButtonCallback, imageh});
% bh2.Position(2) = 80;
% bh3 = uicontrol('Parent',fh,'Style','pushbutton','String','free','Callback',{@freeButtonCallback, imageh});
% bh3.Position(2) = 110;
% 
% % button callback function
% function paintButtonCallback(obj,~,imageh)
% if isempty(obj.Tag)
%     imageh.ButtonDownFcn = @paintMode;
%     obj.Tag = 'on';
% else
%     imageh.ButtonDownFcn = '';
%     obj.Tag = '';
% end
% 
%     function paintMode(~,~)
%         [x,y] = ginput(1);
% 
%         % round the values so they can be used for indexing.
%         x = round(x);
%         y = round(y);
% 
%         % make sure the values do not go outside the image.
%         s = size(imageh.CData);
%         if x > s(2) || y > s(1) || x < 1 || y < 1
%             return
%         end
% 
%         % make the selected pixel white.
%         imageh.CData(y,x) = true;
%     end
% end
% 
% % button callback function
% function lineButtonCallback(~,~,imageh)
% % take two points at a time
% [x,y] = ginput(2);
% 
% % make sure the values do not go outside the image.
% s = size(imageh.CData);
% if any(x > s(2)+0.5 | y > s(1)+0.5 | x < 0.5 | y < 0.5) || (diff(x) == 0 && diff(y) == 0)
%     return
% end
% 
% % find all pixels on the line xy
% ind = findLine(size(imageh.CData),x,y);
% 
% % make the selected pixel white.
% imageh.CData(ind) = true;
% end
% 
% function ind = findLine(s,x,y)
% % Find all pixels that lie between points defined by [x(1),y(1)] and [x(2),y(2)].
% 
% supersampling = 1.2;
% [x,y,~] = improfile(s,round(x),round(y),max([diff(x);diff(y)])*supersampling);
% ind = sub2ind(s,round(x),round(y));
% end
% 
% % button callback function
% function line2ButtonCallback(~,~,imageh)
% % take two points at a time
% h = drawline;
% ind = h.createMask;
% delete(h);
% 
% % make the selected pixel white.
% imageh.CData(ind) = true;
% end
% 
% % button callback function
% function freeButtonCallback(~,~,imageh)
% % take two points at a time
% h = drawfreehand;
% x = h.Position(:,1);
% y = h.Position(:,2);
% delete(h);
% 
% ind = sub2ind(size(imageh.CData),round(y),round(x));
% 
% % make the selected pixel white.
% imageh.CData(ind) = true;
% end

%%%%%

% fontSize = 12;
% 
% % Read in a standard MATLAB gray scale demo image.
% folder = fullfile(matlabroot, '\toolbox\images\imdemos');
% baseFileName = 'cameraman.tif';
% % Get the full filename, with path prepended.
% fullFileName = fullfile(folder, baseFileName);
% % Check if file exists.
% if ~exist(fullFileName, 'file')
%   % File doesn't exist -- didn't find it there.  Check the search path for it.
%   fullFileName = baseFileName; % No path this time.
%   if ~exist(fullFileName, 'file')
%     % Still didn't find it.  Alert user.
%     errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
%     uiwait(warndlg(errorMessage));
%     return;
%   end
% end
% % grayImage = imread(fullFileName);
% grayImage = 0.5*ones(1000,1000);
% imshow(grayImage, []);
% axis on;
% title('Original Grayscale Image', 'FontSize', fontSize);
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
% message = sprintf('Left click and hold to begin drawing a freehand path.\nSimply lift the mouse button to finish.\nDRAW FAST!!!');
% uiwait(msgbox(message));
% % User draws curve on image here.
% hFH = imfreehand();
% % Get the xy coordinates of where they drew.
% xy = hFH.getPosition
% % get rid of imfreehand remnant.
% delete(hFH);
% % Overlay what they drew onto the image.
% hold on; % Keep image, and direction of y axis.
% xCoordinates = xy(:, 1);
% yCoordinates = xy(:, 2);
% plot(xCoordinates, yCoordinates, 'ro', 'LineWidth', 2, 'MarkerSize', 10);
% caption = sprintf('Original Grayscale Image.\nPoints may not lie on adjacent pixels, depends on your speed of drawing!');
% title(caption, 'FontSize', fontSize);
% % Ask user if they want to burn the line into the image.
% promptMessage = sprintf('Do you want to burn the line into the image?');
% titleBarCaption = 'Continue?';
% button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
% if strcmpi(button, 'Yes')
%   cla;
%   hold off;
%   for k = 1 : length(xCoordinates)
%     row = int32(yCoordinates(k));
%     column = int32(xCoordinates(k));
%     
%     row(row<1)=1;
%     column(column<1)=1;
%     
%     grayImage(row, column) = 255;
%   end
%   imshow(grayImage, []);
%   axis on;
%   caption = sprintf('Grayscale Image with Burned In Curve.\nPoints may not lie on adjacent pixels, depends on your speed of drawing!');
%   title(caption, 'FontSize', fontSize);
% end
% % Ask user if they want to interpolate the line to get the "in-between" points that are missed..
% % promptMessage = sprintf('Do you want to interpolate the curve into intervening pixels?');
% % titleBarCaption = 'Continue?';
% % button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
% % if strcmpi(button, 'Cancel')
% %   return;
% % end
% xCoordinates = xy(:, 1);
% yCoordinates = xy(:, 2);
% numberOfKnots = length(xCoordinates);
% % Close gaps that you get when you draw too fast.
% % Use splines to interpolate a smoother curve,
% % with 10 times as many points,
% % that goes exactly through the same data points.
% samplingRateIncrease = 10;
% newXSamplePoints = linspace(1, numberOfKnots, numberOfKnots * samplingRateIncrease);
% % smoothedY = spline(xCoordinates, yCoordinates, newXSamplePoints);
% % Make the 2D array where the top row is the x coordinates and the bottom row is the y coordinates,
% % but with the exception that the left column and right column is a vector that gives the direction of the slope.
% yy = [0, xCoordinates', 0; 1, yCoordinates', 1]
% pp = spline(1:numberOfKnots, yy); % Get interpolant
% smoothedY = ppval(pp, newXSamplePoints); % Get smoothed y values in the "gaps".
% % smoothedY is a 2D array with the x coordinates in the top row and the y coordinates in the bottom row.
% smoothedXCoordinates = smoothedY(1, :)
% smoothedYCoordinates = smoothedY(2, :)
% % Plot smoothedY and show how the line is
% % smooth, and has no sharp bends.
% hold on; % Don't destroy the first curve we plotted.
% hGreenCurve = plot(smoothedXCoordinates, smoothedYCoordinates, '-g');
% title('Spline Interpolation Demo', 'FontSize', 20);
% % But smoothedXCoordinates and smoothedYCoordinates are not in pixel coordinates, they have fractional values.
% % If you want integer pixel values, you have to round.
% intSmoothedXCoordinates = int32(smoothedXCoordinates)
% intSmoothedYCoordinates = int32(smoothedYCoordinates)
% % But now it's possible that some coordinates will be on the same pixel if that's
% % how they rounded according to how they were located to the nearest integer pixel location.
% % So use diff() to remove elements that have the same x and y values.
% diffX = [1, diff(intSmoothedXCoordinates)];
% diffY = [1, diff(intSmoothedYCoordinates)];
% % Find out where both have zero difference from the prior point.
% bothZero = (diffX==0) & (diffY == 0);
% % Remove those from the arrays.
% finalX = intSmoothedXCoordinates(~bothZero);
% finalY = intSmoothedYCoordinates(~bothZero);
% % Now remove the green line.
% delete(hGreenCurve);
% % Plot the final coordinates.
% hGreenCurve = plot(finalX, finalY, '-y');
% h = fill(finalX,finalY,'w');

%%%%%

% img = imread('test.png');
% handles.axes = axes;
% imshow(img, 'Parent', handles.axes);
% axis(handles.axes, 'on');
% % yLimits = get(handles.axes,'YLim');   % If you have the handle already, ...
% % yTicks = yLimits(2)-get(handles.axes,'YTick');
% % set(handles.axes,'YTickLabel',num2str(yTicks.'));
% size(img);
% 
% grayImage = rgb2gray(img);
% figure()
% imshow(grayImage)
% sz = size(grayImage);
% 
% xg = 1:sz(1);
% yg = 1:sz(2);
% F = griddedInterpolant({xg,yg},double(grayImage),'nearest');
% 
% xq = (0:1/10:sz(1))';
% yq = (0:1/10:sz(2))';
% vq = uint8(F({xq,yq}));
% vq = vq(1:end-1,1:end-1);
% figure()
% imshow(vq)
% 
% vq = rescale(vq,1000,12500);

% rescaled = rescale(grayImage,1000,12500);
% rescaled_upsampled = upsample(rescaled,10);
% rescaled_upsampled = rescaled_upsampled.';
% rescaled_upsampled = upsample(rescaled_upsampled,10);
% resample()

% function pixelValues = savePixelValues(rgb)
%     pixelValues = rgb2gray(rgb);
% end

%%%%%

% x = 1;
% 
% tic
% if x
%     disp(x)
% else
%     3
% end
% toc
% 
% tic
% try
%     x
% catch
%     3
% end
% toc

%%%%%

% possibilities = linspace(1,1000,19); % DIFFUSIVITIES
% possibilities = [possibilities possibilities(end)+1];
% test_vals = randsample(possibilities,1000,true);
% test = reshape(test_vals,[100,10]) % DIFFUSIVITIES SAMPLED BY EACH PARTICLE AT EACH TIME POINT
% % test = randi(10,20,10)
% % possibilities = 1:11;
% 
% all_gc = zeros(size(test,1),length(possibilities)-1);
% for i = 1:size(test,1)
%     [GC,GR] = groupcounts(test(i,:)',possibilities,'IncludeEmptyGroups',true);
%     all_gc(i,:) = GC;
% end
% 
% all_gc
% 
% particles = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20'};
% diffs = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S'};
% % tb = table(all_gc)
% 
% % tb = array2table(all_gc) %,'RowNames',particles,'VariableNames',diffs,'DimensionNames',{'Diffusivities','Particles'})
% 
% % heatmap(1:20,possibilities,all_gc)
% 
% h = heatmap(all_gc);
% h.XDisplayLabels = diffs;
% h.YDisplayLabels = particles;

%%%%%

% x = -25:0.1:25;
% mu = 0;
% sigmas = 1:0.1:5;
% all_p = zeros(1,length(x));
% 
% figure()
% for sigma = sigmas
%     p = (1/(sigma*sqrt(2*pi))).*exp(-(1/2).*((x-mu)./sigma).^2);
%     all_p = all_p + p;
%     plot(x,p)
%     hold on
% end
% title('All Gaussian Curves (\mu=0, \sigma=[1:0.1:5])')
% xlabel('x')
% ylabel('y')
% 
% f = fit(x',all_p','gauss1');
% 
% figure()
% plot(f,x,all_p)
% title('Sum of Gaussian Curves and Gaussian Fit')

%%%%%

% msds = [1:5;4:8;5:9];
% dtaus = 1:5;
% 
% errs = std(msds);
% 
% avg_msd = mean(msds);
% 
% plot(dtaus,msds)
% hold on
% errorbar(dtaus,avg_msd,errs)
% 
% hrs = var(msds)./(avg_msd.^2)

%%%%%

% iters = 5;
% len = 1e5;
% tests = zeros(iters,len);
% 
% for i = 1:iters
%     test = normrnd(0,i,1,len);
%     tests(i,:) = test;
%     
% %     figure()
%     hist_object = histogram(test,'Normalization','pdf');
%     fit_var = fitdist(test','Normal');
%     eval_vals = (hist_object.BinEdges(1)-0.2:0.001:hist_object.BinEdges(end)+0.2);
%     fit_pdf = pdf(fit_var,eval_vals);
%     hold on
%     plot(eval_vals,fit_pdf,'LineWidth',2)
%     
%     clearvars test
% end
% 
% idxs = randi(iters,[1 len]);
% tests_sampled = tests(sub2ind(size(tests),idxs,1:len));
% 
% figure()
% hist_object = histogram(tests_sampled,'Normalization','pdf');
% fit_var = fitdist(tests_sampled','Normal');
% eval_vals = (hist_object.BinEdges(1)-0.2:0.001:hist_object.BinEdges(end)+0.2);
% fit_pdf = pdf(fit_var,eval_vals);
% hold on
% plot(eval_vals,fit_pdf,'LineWidth',2)

%%%%%

% % base_time = 0.1;
% % alpha = 0.9;
% % 
% % lower = 0.1;
% % upper = 1;
% iters = 1e3;
% values = zeros(1,iters);
% for i = 1:iters
% %     num = lower + (upper-lower) .* rand;
% %     values(i) = (base_time^alpha)./(rand.^(1+alpha));
%     values(i) = gprnd(0.1,10,0);
% end
% 
% values = values(values<1e4);
% h = histogram(values,'Normalization','pdf');
% % h.BinWidth = 0.1;
% 
% % hold on
% % tau = 0.1:00.1:10;
% % wait_dist = (base_time^alpha)./(tau.^(1+alpha));
% % 
% % % wait_dist = 1./(tau.^(1+alpha));
% % 
% % plot(tau,wait_dist)
% % cumtrapz(tau,wait_dist)
% 
% % area = sum(h.Values)*h.BinWidth

%%%%%

% time_pts = 10000;
% dt_min = 0.1; %seconds
% dt_max = 1e4; %seconds (approx. 2.77 hours)
% log_dt_min = log10(dt_min);
% log_dt_max = log10(dt_max);
% 
% r = log_dt_min + (log_dt_max-log_dt_min).*rand([time_pts,1]);
% times = 10.^r;
% 
% % min(times)
% % max(times)
% % 
% % histogram(times)
% 
% n = 100; %number of particles
% loc_x = zeros(n,time_pts); %particles positions
% 
% for i = 1:n
%     for j = 2:time_pts
%         dx = normrnd(0,1);
%         loc_x(i,j) = loc_x(i,j-1) + times(j)*dx;
%     end
%     disp(i)
% end
% 
% cumulative_times = cumsum(times);
% 
% for i = 1:n
%     plot(cumulative_times,loc_x(i,:))
%     hold on
% end

%%%%%

% num_time_pts = 1e5;
% dim_x = 1e10;
% dim_y = 1e10;
% s = spalloc(dim_y,dim_x,num_time_pts);

%%%%%

viscosity = abs(normrnd(8.9e-4,1e-3,[1 1e5]));
histogram(viscosity)

%%%%%

% tic
% viscosity = zeros(1,1e5);
% for i = 1:1e5
%     viscosity(i) = i;
% end
% toc

%%%%%

% test = betarnd(5,1,[1 1e8]);
% histogram(test)

%%%%%














