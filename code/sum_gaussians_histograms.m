% Combination of Gaussian Curves

clearvars

%%% STORAGE VARIABLE SET-UP %%%
all_displacement_storage_x_1dt = [];
all_displacement_storage_y_1dt = [];

all_displacement_storage_x_10dt = [];
all_displacement_storage_y_10dt = [];

all_displacement_storage_x_100dt = [];
all_displacement_storage_y_100dt = [];

%%% DATA LOADING %%%
diffusivity_list = [2500,3611,4722,5833,6944,8056,9167,10278,11389,12500];
for d = diffusivity_list
    clearvars temp_data_1dt path_str_10dt path_str_100dt
    
    path_str_1dt = strcat(['/Users/AlexBecker/Desktop/Work/University/WagnerLab2021/results/other/sum_gaussians_test/' num2str(d) '/histograms_all_displacements_1dt.mat']);
    temp_data_1dt = load(path_str_1dt);
    all_displacement_storage_x_1dt = [all_displacement_storage_x_1dt temp_data_1dt.all_displacement_storage_x];
    all_displacement_storage_y_1dt = [all_displacement_storage_y_1dt temp_data_1dt.all_displacement_storage_y];
    
    path_str_10dt = strcat(['/Users/AlexBecker/Desktop/Work/University/WagnerLab2021/results/other/sum_gaussians_test/' num2str(d) '/histograms_all_displacements_10dt.mat']);
    temp_data_10dt = load(path_str_10dt);
    all_displacement_storage_x_10dt = [all_displacement_storage_x_10dt temp_data_10dt.all_displacement_storage_x];
    all_displacement_storage_y_10dt = [all_displacement_storage_y_10dt temp_data_10dt.all_displacement_storage_y];
        
    path_str_100dt = strcat(['/Users/AlexBecker/Desktop/Work/University/WagnerLab2021/results/other/sum_gaussians_test/' num2str(d) '/histograms_all_displacements_100dt.mat']);
    temp_data_100dt = load(path_str_100dt);
    all_displacement_storage_x_100dt = [all_displacement_storage_x_100dt temp_data_100dt.all_displacement_storage_x];
    all_displacement_storage_y_100dt = [all_displacement_storage_y_100dt temp_data_100dt.all_displacement_storage_y];
    
end

%%% HISTOGRAMS: 1dt %%%
% x-direction displacements histogram
figure()
hist_obj_x = histogram(all_displacement_storage_x_1dt, 'Normalization', 'pdf');
% fit_data_x = fitdist(all_displacement_storage_x_1dt','Normal'); %obtain fit data
% fit_data_x = fitdist(all_displacement_storage_x_1dt','Logistic'); %obtain fit data
% fit_data_x = fitdist(all_displacement_storage_x_1dt(all_displacement_storage_x_1dt>0)','Exponential'); %obtain fit data
fit_data_x = fitdist(all_displacement_storage_x_1dt(all_displacement_storage_x_1dt>0)','Gamma'); %obtain fit data

eval_vals_x = (hist_obj_x.BinEdges(1)-20:0.1:hist_obj_x.BinEdges(end)+20);
fit_data_pdf_x = pdf(fit_data_x,eval_vals_x); %compute the corresponding PDF
hold on
plot(eval_vals_x,fit_data_pdf_x,'LineWidth',2) %overlay the PDF on top of the histogram

title('Step Size Distribution')
xlabel('\Deltax [10^{-2}\mum]')
ylabel('P(\Deltax, \Delta\tau=1 time point)')

% fit_legend_x = strcat(['Mean = ' num2str(fit_data_x.mu) ', Std. Dev. = ' num2str(fit_data_x.sigma)]);
% legend('Distribution',fit_legend_x)
% fit_legend_x = strcat(['Mean = ' num2str(fit_data_x.mu)]);
% legend('Distribution',fit_legend_x)
fit_legend_x = strcat(['Shape = ' num2str(fit_data_x.a) ', Scale = ' num2str(fit_data_x.b)]);
legend('Distribution',fit_legend_x)

% y-direction displacements histogram
figure()
hist_obj_y = histogram(all_displacement_storage_y_1dt, 'Normalization', 'pdf');
% fit_data_y = fitdist(all_displacement_storage_y_1dt','Normal'); %obtain fit data
% fit_data_y = fitdist(all_displacement_storage_y_1dt','Logistic'); %obtain fit data
% fit_data_y = fitdist(all_displacement_storage_y_1dt(all_displacement_storage_y_1dt>0)','Exponential'); %obtain fit data
fit_data_y = fitdist(all_displacement_storage_y_1dt(all_displacement_storage_y_1dt>0)','Gamma'); %obtain fit data

eval_vals_y = (hist_obj_y.BinEdges(1)-20:0.1:hist_obj_y.BinEdges(end)+20);
fit_data_pdf_y = pdf(fit_data_y,eval_vals_y); %compute the corresponding PDF
hold on
plot(eval_vals_y,fit_data_pdf_y,'LineWidth',2) %overlay the PDF on top of the histogram

title('Step Size Distribution')
xlabel('\Deltay [10^{-2}\mum]')
ylabel('P(\Deltay, \Delta\tau=1 time point)')

% fit_legend_y = strcat(['Mean = ' num2str(fit_data_y.mu) ', Std. Dev. = ' num2str(fit_data_y.sigma)]);
% legend('Distribution',fit_legend_y)
% fit_legend_y = strcat(['Mean = ' num2str(fit_data_y.mu)]);
% legend('Distribution',fit_legend_y)
fit_legend_y = strcat(['Shape = ' num2str(fit_data_y.a) ', Scale = ' num2str(fit_data_y.b)]);
legend('Distribution',fit_legend_y)

%%% HISTOGRAMS: 10dt %%%
% x-direction displacements histogram
figure()
hist_obj_x = histogram(all_displacement_storage_x_10dt, 'Normalization', 'pdf');
% fit_data_x = fitdist(all_displacement_storage_x_10dt','Normal'); %obtain fit data
% fit_data_x = fitdist(all_displacement_storage_x_10dt','Logistic'); %obtain fit data
% fit_data_x = fitdist(all_displacement_storage_x_10dt(all_displacement_storage_x_10dt>0)','Exponential'); %obtain fit data
fit_data_x = fitdist(all_displacement_storage_x_10dt(all_displacement_storage_x_10dt>0)','Gamma'); %obtain fit data

eval_vals_x = (hist_obj_x.BinEdges(1)-20:0.1:hist_obj_x.BinEdges(end)+20);
fit_data_pdf_x = pdf(fit_data_x,eval_vals_x); %compute the corresponding PDF
hold on
plot(eval_vals_x,fit_data_pdf_x,'LineWidth',2) %overlay the PDF on top of the histogram

title('Step Size Distribution')
xlabel('\Deltax [10^{-2}\mum]')
ylabel('P(\Deltax, \Delta\tau=10 time points)')

% fit_legend_x = strcat(['Mean = ' num2str(fit_data_x.mu) ', Std. Dev. = ' num2str(fit_data_x.sigma)]);
% legend('Distribution',fit_legend_x)
% fit_legend_x = strcat(['Mean = ' num2str(fit_data_x.mu)]);
% legend('Distribution',fit_legend_x)
fit_legend_x = strcat(['Shape = ' num2str(fit_data_x.a) ', Scale = ' num2str(fit_data_x.b)]);
legend('Distribution',fit_legend_x)

% y-direction displacements histogram
figure()
hist_obj_y = histogram(all_displacement_storage_y_10dt, 'Normalization', 'pdf');
% fit_data_y = fitdist(all_displacement_storage_y_10dt','Normal'); %obtain fit data
% fit_data_y = fitdist(all_displacement_storage_y_10dt','Logistic'); %obtain fit data
% fit_data_y = fitdist(all_displacement_storage_y_10dt(all_displacement_storage_y_10dt>0)','Exponential'); %obtain fit data
fit_data_y = fitdist(all_displacement_storage_y_10dt(all_displacement_storage_y_10dt>0)','Gamma'); %obtain fit data

eval_vals_y = (hist_obj_y.BinEdges(1)-20:0.1:hist_obj_y.BinEdges(end)+20);
fit_data_pdf_y = pdf(fit_data_y,eval_vals_y); %compute the corresponding PDF
hold on
plot(eval_vals_y,fit_data_pdf_y,'LineWidth',2) %overlay the PDF on top of the histogram

title('Step Size Distribution')
xlabel('\Deltay [10^{-2}\mum]')
ylabel('P(\Deltay, \Delta\tau=10 time points)')

% fit_legend_y = strcat(['Mean = ' num2str(fit_data_y.mu) ', Std. Dev. = ' num2str(fit_data_y.sigma)]);
% legend('Distribution',fit_legend_y)
% fit_legend_y = strcat(['Mean = ' num2str(fit_data_y.mu)]);
% legend('Distribution',fit_legend_y)
fit_legend_y = strcat(['Shape = ' num2str(fit_data_y.a) ', Scale = ' num2str(fit_data_y.b)]);
legend('Distribution',fit_legend_y)

%%% HISTOGRAMS: 100dt %%%
% x-direction displacements histogram
figure()
hist_obj_x = histogram(all_displacement_storage_x_100dt, 'Normalization', 'pdf');
% fit_data_x = fitdist(all_displacement_storage_x_100dt','Normal'); %obtain fit data
% fit_data_x = fitdist(all_displacement_storage_x_100dt','Logistic'); %obtain fit data
% fit_data_x = fitdist(all_displacement_storage_x_100dt(all_displacement_storage_x_100dt>0)','Exponential'); %obtain fit data
fit_data_x = fitdist(all_displacement_storage_x_100dt(all_displacement_storage_x_100dt>0)','Gamma'); %obtain fit data

eval_vals_x = (hist_obj_x.BinEdges(1)-20:0.1:hist_obj_x.BinEdges(end)+20);
fit_data_pdf_x = pdf(fit_data_x,eval_vals_x); %compute the corresponding PDF
hold on
plot(eval_vals_x,fit_data_pdf_x,'LineWidth',2) %overlay the PDF on top of the histogram

title('Step Size Distribution')
xlabel('\Deltax [10^{-2}\mum]')
ylabel('P(\Deltax, \Delta\tau=100 time points)')

% fit_legend_x = strcat(['Mean = ' num2str(fit_data_x.mu) ', Std. Dev. = ' num2str(fit_data_x.sigma)]);
% legend('Distribution',fit_legend_x)
% fit_legend_x = strcat(['Mean = ' num2str(fit_data_x.mu)]);
% legend('Distribution',fit_legend_x)
fit_legend_x = strcat(['Shape = ' num2str(fit_data_x.a) ', Scale = ' num2str(fit_data_x.b)]);
legend('Distribution',fit_legend_x)

% y-direction displacements histogram
figure()
hist_obj_y = histogram(all_displacement_storage_y_100dt, 'Normalization', 'pdf');
% fit_data_y = fitdist(all_displacement_storage_y_100dt','Normal'); %obtain fit data
% fit_data_y = fitdist(all_displacement_storage_y_100dt','Logistic'); %obtain fit data
% fit_data_y = fitdist(all_displacement_storage_y_100dt(all_displacement_storage_y_100dt>0)','Exponential'); %obtain fit data
fit_data_y = fitdist(all_displacement_storage_y_100dt(all_displacement_storage_y_100dt>0)','Gamma'); %obtain fit data

eval_vals_y = (hist_obj_y.BinEdges(1)-20:0.1:hist_obj_y.BinEdges(end)+20);
fit_data_pdf_y = pdf(fit_data_y,eval_vals_y); %compute the corresponding PDF
hold on
plot(eval_vals_y,fit_data_pdf_y,'LineWidth',2) %overlay the PDF on top of the histogram

title('Step Size Distribution')
xlabel('\Deltay [10^{-2}\mum]')
ylabel('P(\Deltay, \Delta\tau=100 time points)')

% fit_legend_y = strcat(['Mean = ' num2str(fit_data_y.mu) ', Std. Dev. = ' num2str(fit_data_y.sigma)]);
% legend('Distribution',fit_legend_y)
% fit_legend_y = strcat(['Mean = ' num2str(fit_data_y.mu)]);
% legend('Distribution',fit_legend_y)
fit_legend_y = strcat(['Shape = ' num2str(fit_data_y.a) ', Scale = ' num2str(fit_data_y.b)]);
legend('Distribution',fit_legend_y)