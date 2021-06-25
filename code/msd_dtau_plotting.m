function msd_dtau_plotting(n,time_pts,data_matrix,boundary_collision,conversion_factor,msd_dtau_log)
% MSD_DTAU_PLOTTING Plot the time-averaged MSD(delta-tau) curve for each
% particle as well as the ensemble-averaged MSD(delta-tau) curve.

%%% PLOT SET-UP %%%
time_fraction = 1/3; %fraction of the total number of time points over which to generate the curve
delta_taus = 1:round(time_pts*time_fraction); %time point intervals for displacement measurements
sqd_dispmnts_lag_time = zeros(n,time_pts,size(delta_taus,2)); %storage for displacements at each time point interval
rescaled_time = 0.1:0.1:time_pts*conversion_factor; %time in seconds used for more realistic plotting

% Compute the displacements for the given delta-tau values (multiples of delta-t)
counter_msd_tau = 1;
for dt = delta_taus
    for i = 1:n
        for j = 1:time_pts-dt
            tamsd_plot_displacement_x = data_matrix(i,j+dt,1) - data_matrix(i,j,1);
            tamsd_plot_displacement_y = data_matrix(i,j+dt,2) - data_matrix(i,j,2);
            sqd_dispmnts_lag_time(i,j,counter_msd_tau) = tamsd_plot_displacement_x^2 + tamsd_plot_displacement_y^2;
        end
    end
    counter_msd_tau = counter_msd_tau + 1;
end

% Remove erroneous displacements that result from a particle striking the boundary
for i = 1:n
    for j = 1:size(delta_taus,2)
        if boundary_collision(i) == 1 %only modify the displacement data if the given particle strikes the boundary
            start_of_trailing_zeros = find(sqd_dispmnts_lag_time(i,:,j),1,'last') + 1;
            try
                % Remove the appropriate number of erroneous displacements
                sqd_dispmnts_lag_time(i,(start_of_trailing_zeros-delta_taus(j)):start_of_trailing_zeros-1,j) = 0; 
            catch
                % If the particle stopped so early that all displacements are erroneous, set the enter set of displacements to zero
                sqd_dispmnts_lag_time(i,:,j) = 0;
            end
        end
    end
end

%%% PLOTTING %%%
% Set up the plot for the time-averaged squared displacements versus lag time curves
sdlt_plotting = zeros(n,size(delta_taus,2));
for i = 1:n
    for j = 1:size(delta_taus,2)
        current_particle_tau = sqd_dispmnts_lag_time(i,:,j); %isolate one particle
        start_of_trailing_zeros = find(sqd_dispmnts_lag_time(i,:,j),1,'last') + 1; %find the index of the two consecutive zeros in sqd_dispmnts_lag_time for the given multiple of delta-t
        current_particle_tau(start_of_trailing_zeros:end) = []; %remove all trailing zeros from the data matrix
        sdlt_plotting(i,j) = mean(current_particle_tau); %compute the mean for the given particle at the given time lag
    end
end

% Plot the time-averaged squared displacements versus lag time curves
figure()
for i = 1:n
    clearvars sdlt_plotting_particle
    sdlt_plotting_particle = sdlt_plotting(i,:);
    sdlt_plotting_particle = no_trailing_zeros(sdlt_plotting_particle); %remove trailing (excess) zeros if a given walk ended early
    sdlt_plotting_particle = 10^-4.*(sdlt_plotting_particle/conversion_factor);
    plot(rescaled_time(1:length(sdlt_plotting_particle)),sdlt_plotting_particle)
    hold on
    
    % Line of best fit for each individual particle's curve
    [Dfit, alphafit] = msd_dtau_fitting(0.01,1,delta_taus(1:length(sdlt_plotting_particle)),sdlt_plotting_particle);
    diffusivity_disp_str = strcat(['Diffusion constant from the line of best fit (particle #' num2str(i) '): ' num2str(Dfit)]);
    alpha_disp_str = strcat(['Alpha from the line of best fit (particle #' num2str(i) '): ' num2str(alphafit)]);
    disp(diffusivity_disp_str)
    disp(alpha_disp_str)
end

% Convert NaN values (resulting from entirely-zero walks) to zeros to avoid breaking the MSD(delta-tau) curve
sdlt_plotting(isnan(sdlt_plotting)) = 0;

% Set up the plot for the corresponding time-averaged MSD(delta-tau) curve
msd_tau_plotting = zeros(1,size(delta_taus,2));
for i = 1:size(delta_taus,2)
    sdlt_plotting_no_zeros = sdlt_plotting(:,i);
    sdlt_plotting_no_zeros(sdlt_plotting_no_zeros==0) = [];
    msd_tau_plotting(i) = mean(sdlt_plotting_no_zeros);
end

% Plot the corresponding time-averaged MSD(delta-tau) curve
msd_tau_plotting = 10^-4.*(msd_tau_plotting/conversion_factor);
plot(rescaled_time(1:length(msd_tau_plotting)),msd_tau_plotting,'LineWidth',3,'Color','red')

% Plot title and axis labels
title('Time-Averaged Squared Displacement vs. Lag Time')
xlabel('Lag Time \Delta\tau [s]')
ylabel('Time-Averaged MSD(\Delta\tau) [\mum^2]')

% Scale the axes logarithmically if requested
if msd_dtau_log == 1
    set(gca,'XScale','log');
    set(gca,'YScale','log');
end

%%% LINE OF BEST FIT %%%
[Dfit, alphafit] = msd_dtau_fitting(0.01,1,delta_taus,msd_tau_plotting);
diffusivity_disp_str = strcat(['Diffusion constant from the line of best fit: ' num2str(Dfit)]);
alpha_disp_str = strcat(['Alpha from the line of best fit: ' num2str(alphafit)]);
plot((delta_taus*conversion_factor),4.*Dfit.*delta_taus.^alphafit,'Linewidth',3,'Color','b','LineStyle',':')
disp(diffusivity_disp_str)
disp(alpha_disp_str)

%%% SEGMENTED LINES OF BEST FIT %%%
log_times = log10(rescaled_time(1:length(sdlt_plotting_particle)));
segments = rescaled_time(log_times==round(log_times)); %get all the powers of ten seconds

% For each lag time segment bounded by a power of ten seconds, compute the parameters of a line of best fit
for i = 1:length(segments)
    try
        [Dfit, alphafit] = msd_dtau_fitting(0.01,1,delta_taus(find(rescaled_time==segments(i)):find(rescaled_time==segments(i+1))),msd_tau_plotting(find(rescaled_time==segments(i)):find(rescaled_time==segments(i+1))));
        diffusivity_disp_str = strcat(['Diffusion constant from the line of best fit (' num2str(segments(i)) 's - ' num2str(segments(i+1)) 's' '): ' num2str(Dfit)]);
        alpha_disp_str = strcat(['Alpha from the line of best fit: (' num2str(segments(i)) 's - ' num2str(segments(i+1)) 's' '): ' num2str(alphafit)]);
        plot((delta_taus(find(rescaled_time==segments(i)):find(rescaled_time==segments(i+1)))*conversion_factor),4.*Dfit.*delta_taus(find(rescaled_time==segments(i)):find(rescaled_time==segments(i+1))).^alphafit,'Linewidth',3,'Color','k','LineStyle','--')
    catch %used for the final segment (from the final segment value in "segments" to the end of the lag time range)
        if length(log_times)-(segments(end)/conversion_factor) < 10 %ignore this final segment if there are too few data points for it to be representative of a trend (the minimum is hard-coded to ten data points)
            disp('NOTE. The final bracket is too small to attempt a slope and intercept estimate.')
            break
        else
            [Dfit, alphafit] = msd_dtau_fitting(0.01,1,delta_taus(find(rescaled_time==segments(i)):end),msd_tau_plotting(find(rescaled_time==segments(i)):end));
            diffusivity_disp_str = strcat(['Diffusion constant from the line of best fit (' num2str(segments(i)) 's - end): ' num2str(Dfit)]);
            alpha_disp_str = strcat(['Alpha from the line of best fit: (' num2str(segments(i)) 's - end): ' num2str(alphafit)]);
            plot((delta_taus(find(rescaled_time==segments(i)):end)*conversion_factor),4.*Dfit.*delta_taus(find(rescaled_time==segments(i)):end).^alphafit,'Linewidth',3,'Color','k','LineStyle','--')
        end
    end
    
    disp(diffusivity_disp_str)
    disp(alpha_disp_str)
    
end

end