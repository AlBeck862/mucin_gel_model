function walk_visualization(n,lattice,data_matrix)
% WALK_VISUALIZATION Plot the random walk of each simulated particle.

if length(unique(lattice)) == 1
    lattice_visualization = rescale(lattice,1,1);
else
    lattice_visualization = rescale(lattice,0,1);
end

lattice_visualization = fliplr(rot90(lattice_visualization,-1));

for i = 1:n
    pause(3)
    
    figure()
    imshow(lattice_visualization)
    hold on
    
    % Remove excess zeros (in the event the particle struck the boundary)
    current_particle_x = data_matrix(i,:,1);
    current_particle_x(current_particle_x==0) = [];
    
    % Remove excess zeros (in the event the particle struck the boundary)
    current_particle_y = data_matrix(i,:,2);
    current_particle_y(current_particle_y==0) = [];
    
    plot(current_particle_x,current_particle_y)                                                         %plot the entire trajectory
    plot(current_particle_x(1),current_particle_y(1),'>g','MarkerFaceColor','g','MarkerSize',10)        %mark the start point
    plot(current_particle_x(end),current_particle_y(end),'sr','MarkerFaceColor','r','MarkerSize',10)    %mark the end point
    legend('Particle Trajectory','Start','End')
    title_str = strcat(['Particle Trajectory (Particle #' num2str(i) ')']);
    title(title_str)
    
    file_str = strcat(['/temp_results/walks/walk_particle' num2str(i) '.png']);
    saveas(gcf,[pwd file_str]);
    
    clearvars current_particle_x current_particle_y %clear variables for the next loop iteration
end

end