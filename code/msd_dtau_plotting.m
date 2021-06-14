function msd_dtau_plotting(n,time_pts,data_matrix,boundary_collision)
% MSD_DTAU_PLOTTING Plot the time-averaged MSD(delta-tau) curve for each
% particle as well as the ensemble-averaged MSD(delta-tau) curve.

% Time-averaged MSD(delta-tau) plot
delta_taus = 1:(time_pts/10); %time point intervals for displacement measurements (given by the first 10% of time points)
sqd_dispmnts_lag_time = zeros(n,time_pts,size(delta_taus,2)); %storage for displacements at each time point interval

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
    plot(sdlt_plotting(i,:))
    hold on
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
plot(msd_tau_plotting,'LineWidth',2,'Color','red')
title('Time-Averaged Squared Displacement vs. Lag Time')
xlabel('Lag Time \Delta\tau [simulation time points]')
ylabel('Time-Averaged MSD(\Delta\tau) [(10^{-2}\mum)^2]')

% Fetch the line of best fit for the MSD(tau)
msd_tau_fit = polyfit(log10(delta_taus(1:end/2)),log10(msd_tau_plotting(1:end/2)),1);
disp(msd_tau_fit)

end