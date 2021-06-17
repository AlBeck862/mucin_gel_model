function stepsize_dt_plotting(n,data_matrix,moving_avg_kernel)
% STEPSIZE_IN_TIME Plot the step size over time for each particle in both
% the x-direction and the y-direction.

for i = 1:n
    current_particle_x = no_trailing_zeros(data_matrix(i,:,1));
    current_particle_y = no_trailing_zeros(data_matrix(i,:,2));
    
    stepsizes_x = abs(diff(current_particle_x));
    stepsizes_y = abs(diff(current_particle_y));
    
    figure()
    plot(movmean(stepsizes_x,moving_avg_kernel))
    title('Step-Size (x-direction) vs. Absolute Time')
    xlabel('Absolute Time [time points]')
    ylabel('Step-Size |x| [(10^{-2}\mum)]')
    
    figure()
    plot(movmean(stepsizes_y,moving_avg_kernel))
    title('Step-Size (y-direction) vs. Absolute Time')
    xlabel('Absolute Time [time points]')
    ylabel('Step-Size |y| [(10^{-2}\mum)]')
    
end

end