function start_coords = start_location(random_start,lattice_x,lattice_y)
% START_COORDS (explanation)

if random_start == 1
    % A margin must be set to avoid starting too close to the boundary
    margin_percent_x = 10; %the particle won't start within 10% of the lattice's x-length from the boundary on both sides (e.g.: lattice_x=10000 --> the particle's starting x-coordinate can be between 1000 and 9000)
    start_margin_x = round(lattice_x/margin_percent_x);
    min_start_x = start_margin_x;
    max_start_x = lattice_x-start_margin_x;
    
    % A margin must be set to avoid starting too close to the boundary
    margin_percent_y = 10; %the particle won't start within 10% of the lattice's y-length from the boundary on both sides (e.g.: lattice_y=10000 --> the particle's starting y-coordinate can be between 1000 and 9000)
    start_margin_y = round(lattice_y/margin_percent_y);
    min_start_y = start_margin_y;
    max_start_y = lattice_y-start_margin_y;
    
    % Generate a random start location
    % Select a random start location at a sufficient distance from the boundary
    x_start = round(min_start_x + (max_start_x-min_start_x).*rand());
    y_start = round(min_start_y + (max_start_y-min_start_y).*rand());
elseif random_start == -1 %for testing purposes
    %  The particle will start at a specific, hard-coded location
	x_start = 200;
    y_start = 200;
else
    % The particle will start at the center of the lattice
    x_start = round(lattice_x/2);
    y_start = round(lattice_y/2);
    
end

start_coords = [x_start,y_start];
    
end