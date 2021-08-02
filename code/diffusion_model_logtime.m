% Log-time diffusion model
tic

fid = fopen('temp_results/data.txt','wt'); %open a text file to store simulation data

% Define the smallest and largest possible time step sizes
dt_min = 1; %seconds
dt_max = 1e10; %seconds
log_dt_min = log10(dt_min);
log_dt_max = log10(dt_max);

% Basic simulation parameters
num_particles = 5e3; %number of particles to simulate
num_time_pts = 5e4; %number of time points (each separated by a randomly-sized time step)
conversion_factor = 1; %seconds per time point

% Source: Khan and Mason (2014)
% Define the time steps
tic
q = log_dt_min + (log_dt_max-log_dt_min)*rand([num_time_pts,1]);
time_pts = round(10.^q);
time_pts = sort(unique(time_pts));
num_time_pts = length(time_pts); %redefine the number of time steps (in case multiple identical time steps were drawn randomly)
bench_time = toc;
bench_time_disp = strcat(['Time steps defined in ' num2str(bench_time) ' seconds.']);
disp(bench_time_disp)
fprintf(fid,'Time steps defined in %f seconds.\n',bench_time);

% Source: Khan and Mason (2014)
% Generate basic normal distributions
tic
base_dist_x = normrnd(0,1,[num_particles,num_time_pts]);
base_dist_y = normrnd(0,1,[num_particles,num_time_pts]);
bench_time = toc;
bench_time_disp = strcat(['Base distributions generated in ' num2str(bench_time) ' seconds.']);
disp(bench_time_disp)
fprintf(fid,'Base distributions generated in %f seconds.\n',bench_time);

% Define the diffusion coefficient
radius = 1e-6; %m
% viscosity = 8.9e-4; %Pa*s [kg/(m*s)]
viscosity = abs(normrnd(8.9e-4,1e-4,[1 num_time_pts])); %Pa*s [kg/(m*s)]
% viscosity = abs(normrnd(8.9e-4,1e-3,[1 num_time_pts])); %Pa*s [kg/(m*s)]

% viscosity = zeros(1,num_time_pts);
% num_sections = 1e4;
% start_idx = 1;
% end_idx = round(num_time_pts/num_sections);
% for i = 1:num_sections
%     if i == num_sections
%         viscosity(start_idx:end) = 8.9e-4*i;
%         break
%     end
%     
%     viscosity(start_idx:end_idx) = 8.9e-4*i;
%     start_idx = end_idx;
%     end_idx = end_idx + round(num_time_pts/num_sections);
% end

% for i = 1:num_time_pts
%     select = 
%     
% end

kb = 1.38064852e-23; %m^2*kg/(s^2*K)
T = 310; %K (37C)
diffusivity = 1e12.*((kb*T)./(6*pi*viscosity*radius)); %um^2/s
% diffusivity = 1; %um^2/s TEMPORARY OVERRIDE
fprintf(fid,'Particle radius: %f, medium viscosity: %f, diffusivity: %f\n\n',radius,viscosity,diffusivity);

% Source: Khan and Mason (2014)
% Scale the normal distributions according to each time step
tic
mult_factor = sqrt(2*diffusivity.*time_pts'); %multiply by the conversion factor? *****
scaled_dist_x = mult_factor.*base_dist_x;
scaled_dist_y = mult_factor.*base_dist_y;
bench_time = toc;
bench_time_disp = strcat(['Scaled distributions generated in ' num2str(bench_time) ' seconds.']);
disp(bench_time_disp)
fprintf(fid,'Scaled distributions generated in %f seconds.\n',bench_time);

% Compute the overall trajectory of each particle
tic
trajectory_x = cumsum(scaled_dist_x,2);
trajectory_y = cumsum(scaled_dist_y,2);
bench_time = toc;
bench_time_disp = strcat(['Overall trajectories computed in ' num2str(bench_time) ' seconds.']);
disp(bench_time_disp)
fprintf(fid,'Overall trajectories computed in %f seconds.\n',bench_time);

%%% WALK VISUALIZATIONS %%%
% Plot certain randomly-selected trajectories
walk_plot_idxs = sort(unique(randi(num_particles,1,10))); %arbitrarily select up to 10 trajectories to plot
for i = walk_plot_idxs
    figure()
    plot(trajectory_x(i,:),trajectory_y(i,:))
    hold on
    plot(trajectory_x(i,1),trajectory_y(i,1),'>g','MarkerFaceColor','g','MarkerSize',10)        %mark the start point
    plot(trajectory_x(i,end),trajectory_y(i,end),'sr','MarkerFaceColor','r','MarkerSize',10)    %mark the end point
    
    title_str = strcat(['Particle Trajectory (Particle #' num2str(i) ')']);
    title(title_str)
    xlabel('[\mum]')
    ylabel('[\mum]')
    legend('Particle Trajectory','Start','End')
    
    aspect_ratio = get(gca,'DataAspectRatio');
    if aspect_ratio(1) > aspect_ratio(2)
        aspect_ratio(2) = aspect_ratio(1);
    else
        aspect_ratio(1) = aspect_ratio(2);
    end
    set(gca,'DataAspectRatio',aspect_ratio)
end

%%% MSD(dtau) CURVE %%%
% Compute MSD(dtau) for all steps for each particle
% Note: given the logarithmic time step spacing, the dtau value automatically gets larger throughout the simulation.
tic
msd_dtau_all = scaled_dist_x.^2 + scaled_dist_y.^2;

% Compute the overall MSD(dtau)
msd_dtau = mean(msd_dtau_all,1);

% Plot the MSD(dtau) curve on logarithmic axes
figure()
loglog(time_pts,msd_dtau,'Linewidth',3,'Color','k')
title('MSD(\Delta\tau)')
xlabel('\Delta\tau [s]')
ylabel('MSD [\mum^2/s]')

% Fit the MSD(dtau) curve
[Dfit, alphafit] = msd_dtau_fitting(0.01,1,time_pts,msd_dtau);
diffusivity_disp_str = strcat(['Diffusion constant from the line of best fit: ' num2str(Dfit)]);
alpha_disp_str = strcat(['Alpha from the line of best fit: ' num2str(alphafit)]);        
disp(diffusivity_disp_str)
disp(alpha_disp_str)

% Plot the line of best fit
hold on
plot((time_pts*conversion_factor),4.*Dfit.*time_pts.^alphafit,'Linewidth',2,'Color','r','LineStyle',':')
legend('MSD(\Delta\tau)','Line of Best Fit')

bench_time = toc;
bench_time_disp = strcat(['MSD computed in ' num2str(bench_time) ' seconds.']);
disp(bench_time_disp)
fprintf(fid,'MSD computed in %f seconds.\n',bench_time);

%%% HISTOGRAMS %%%
hist_plot_idxs = sort(unique(randi(num_time_pts,1,10))); %arbitrarily select up to 10 time steps to plot
for i = hist_plot_idxs
    temp_data_x = scaled_dist_x(:,i);
    temp_data_y = scaled_dist_y(:,i);
    
    figure()
    hist_obj_x = histogram(temp_data_x,'Normalization','pdf');
    fit_data_x = fitdist(temp_data_x,'Normal');
    eval_vals_x = (hist_obj_x.BinEdges(1):0.01:hist_obj_x.BinEdges(end));
    fit_data_pdf_x = pdf(fit_data_x,eval_vals_x);
    hold on
    plot(eval_vals_x,fit_data_pdf_x,'LineWidth',2)
    
    title('Step Size Distribution')
    xlabel('\Deltax [\mum]')
    hist_y_label_str = strcat(['P(\Deltax, \Delta\tau=' num2str(time_pts(i)) ' seconds)']);
    ylabel(hist_y_label_str)
	fit_legend_x = strcat(['Mean = ' num2str(fit_data_x.mu) ', Std. Dev. = ' num2str(fit_data_x.sigma)]);
    legend('Distribution',fit_legend_x)
    
    figure()
    hist_obj_y = histogram(temp_data_y,'Normalization','pdf');
    fit_data_y = fitdist(temp_data_y,'Normal');
    eval_vals_y = (hist_obj_y.BinEdges(1):0.01:hist_obj_y.BinEdges(end));
    fit_data_pdf_y = pdf(fit_data_y,eval_vals_y);
    hold on
    plot(eval_vals_y,fit_data_pdf_y,'LineWidth',2)
    
    title('Step Size Distribution')
    xlabel('\Deltay [\mum]')
    hist_y_label_str = strcat(['P(\Deltay, \Delta\tau=' num2str(time_pts(i)) ' seconds)']);
    ylabel(hist_y_label_str)
	fit_legend_y = strcat(['Mean = ' num2str(fit_data_y.mu) ', Std. Dev. = ' num2str(fit_data_y.sigma)]);
    legend('Distribution',fit_legend_y)
    
    clearvars temp_data_x temp_data_y
end

bench_time = toc;
bench_time_disp = strcat(['Simulation complete in ' num2str(bench_time) ' seconds.']);
disp(bench_time_disp)
fprintf(fid,'Simulation complete in %f seconds.\n\n',bench_time);

fprintf(fid,'%s\n%s',diffusivity_disp_str,alpha_disp_str);

fclose(fid);