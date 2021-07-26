% Log-time diffusion model

dt_min = 1; %seconds
dt_max = 1e6; %seconds
log_dt_min = log10(dt_min);
log_dt_max = log10(dt_max);

conversion_factor = 1; %seconds per time point
n = 1e4; %number of particles

radius = 1e-6; %m
viscosity = 8.9e-4; %Pa*s [kg/(m*s)]
kb = 1.38064852e-23; %m^2*kg/(s^2*K)
T = 310; %K (37C)
diffusivity = 1e12*(kb*T/(6*pi*viscosity*radius)); %um^2/s
diffusivity = 1; %um^2/s TEMPORARY OVERRIDE

%%% NEW %%%
% Source: Khan and Mason (2014)
num_time_pts = 1e5;
q = log_dt_min + (log_dt_max-log_dt_min)*rand([num_time_pts,1]);
time_pts = round(10.^q);
time_pts = sort(unique(time_pts));
num_time_pts = length(time_pts);

base_dist_x = normrnd(0,1,[n,num_time_pts]);
base_dist_y = normrnd(0,1,[n,num_time_pts]);

mult_factor = sqrt(2*diffusivity*time_pts');
scaled_dist_x = mult_factor.*base_dist_x;
scaled_dist_y = mult_factor.*base_dist_y;

trajectory_x = cumsum(scaled_dist_x,2);
trajectory_y = cumsum(scaled_dist_y,2);

% for i = 1:n
%     figure()
%     plot(trajectory_x(i,:),trajectory_y(i,:))
%     hold on
%     plot(trajectory_x(i,1),trajectory_y(i,1),'>g','MarkerFaceColor','g','MarkerSize',10)        %mark the start point
%     plot(trajectory_x(i,end),trajectory_y(i,end),'sr','MarkerFaceColor','r','MarkerSize',10)    %mark the end point
% end

msd_dtau_all = scaled_dist_x.^2 + scaled_dist_y.^2;
% size(msd_dt_all)

msd_dtau = mean(msd_dtau_all,1);

figure()
loglog(time_pts,msd_dtau)


% HISTOGRAMS: these would be for all particles, display a histogram of all
% displacements for each delta-time (one histogram per column in the
% matrix, aka one histogram per delta-time)

bleeeeh

%%% NEW %%%



%%% OLD STUFF, POSSIBLY TO DISCARD OR RECYCLE %%%
dtau_min = 1; %seconds
dtau_max = 1e6; %seconds
log_dtau_min = log10(dtau_min);
log_dtau_max = log10(dtau_max);

% Source: Khan and Mason (2014)
num_dtaus = 1e2; %number of dtau values for which to compute MSD(dtau) values
r = log_dtau_min + (log_dtau_max-log_dtau_min)*rand([num_dtaus,1]);
dtaus = round(10.^r);
dtaus = sort(unique(dtaus));
num_dtaus = length(dtaus);

loc_x = zeros(n,num_time_pts); %particles x positions
loc_y = zeros(n,num_time_pts); %particles y positions
sigma = sqrt(2*(diffusivity)*conversion_factor);
for i = 1:n
    for j = 2:num_time_pts
        loc_x(i,j) = loc_x(i,j-1) + normrnd(0,sigma);
        loc_y(i,j) = loc_y(i,j-1) + normrnd(0,sigma);
    end
    disp(i)
end

msds_per_particle = zeros(n,num_dtaus);
iteration = 1;
for dtau = dtaus'
    if dtau == 1
        for i = 1:n
             msd_temp = (loc_x(i,2:dtau:end) - loc_x(i,1:dtau:(end-1))).^2 + (loc_y(i,2:dtau:end) - loc_y(i,1:dtau:(end-1))).^2;
             msds_per_particle(i,iteration) = mean(msd_temp);
        end
    else
        for i = 1:n
            msd_temp = (loc_x(i,(dtau+1):dtau:end) - loc_x(i,1:dtau:(end-dtau))).^2 + (loc_y(i,(dtau+1):dtau:end) - loc_y(i,1:dtau:(end-dtau))).^2;
            msds_per_particle(i,iteration) = mean(msd_temp);
        end
    end
    iteration = iteration + 1;
end

msds = mean(msds_per_particle,1);

figure()
plot(dtaus,msds)

figure()
loglog(dtaus,msds)

[Dfit, alphafit] = msd_dtau_fitting(0.01,1,dtaus,msds)