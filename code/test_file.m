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

for i = 1:n
    for j = 1:size(delta_taus,2)
        if boundary_collision(i) == 1 %only modify the displacement data if the given  particle strikes the boundary
            first_zero_idx = find(((sqd_dispmnts_lag_time(i,:,j)==0)+([diff(sqd_dispmnts_lag_time(i,:,j)) 0]==0))==2,1); %find the index of the two consecutive zeros in sqd_dispmnts_lag_time for the given multiple of delta-t
            try
                sqd_dispmnts_lag_time(i,(first_zero_idx-delta_taus(j)):first_zero_idx-1,j) = 0; %remove the appropriate number of erroneous displacements
            catch
                sqd_dispmnts_lag_time(i,:,j) = 0; %remove all displacements if the particle is immobilized at a time point lesser than the value of the time lag
            end
        end
    end
end