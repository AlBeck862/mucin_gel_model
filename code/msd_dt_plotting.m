function msd_dt_plotting(n,time_pts,data_matrix,conversion_factor,multiplier)
% MSD_DT_PLOTTING Plot the MSD(delta-time) curve for each particle as well
% as the average MSD(delta-time) curve.

% Plot set-up
squared_displacements_abs_time = zeros(n,time_pts);
rescaled_time = conversion_factor:conversion_factor:time_pts*conversion_factor; %time in seconds used for more realistic plotting

% Compute the displacements from the start point to each point in the simulation
for i = 1:n
    for j = 1:time_pts
        displacement_x = data_matrix(i,j,1) - data_matrix(i,1,1);
        displacement_y = data_matrix(i,j,2) - data_matrix(i,1,2);
        squared_displacements_abs_time(i,j) = displacement_x^2 + displacement_y^2;
    end
end

% Generate the squared displacement versus absolute time curves
figure()
for i = 1:n
    plot_sdat = squared_displacements_abs_time(i,:);
    plot_sdat(plot_sdat==plot_sdat(end)) = []; %remove the artifically inflated final values (applicable only when a particle was immobilized by the boundary)
    plot_sdat = (1/multiplier).*(plot_sdat./conversion_factor);
    plot(rescaled_time(1:length(plot_sdat)),plot_sdat,'Color','black')
    hold on
    
    clearvars plot_sdat
end

% Generate the corresponding MSD(delta-t) curve
msd_dt_data = zeros(size(squared_displacements_abs_time));
for i = 1:n
    plot_msd_dt = squared_displacements_abs_time(i,:);
    plot_msd_dt(plot_msd_dt==plot_msd_dt(end)) = 0; %replace the artifically inflated final values with zero (applicable only when a particle was immobilized by the boundary)
    msd_dt_data(i,:) = plot_msd_dt;

    clearvars plot_msd_dt
end

% Plot the corresponding MSD(delta-t) curve
msd_dt = sum(msd_dt_data,1)./sum(msd_dt_data~=0,1);
msd_dt = (1/multiplier).*(msd_dt./conversion_factor);
plot(rescaled_time(1:length(msd_dt)),msd_dt,'LineWidth',3,'Color','red')
title('Squared Displacement vs. Absolute Time')
xlabel('Absolute Time \Deltat [seconds]')
ylabel('MSD(\Deltat) [\mum^2]')

saveas(gcf,[pwd '/temp_results/msd_dt.jpeg']);

end