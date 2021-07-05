function cdf = gen_PDF(D,tau,x)
% GEN_PDF Generate a PDF and CDF for a given diffusivity (D), lag time
% (tau), and range of x values (x).

% Define the PDF given the diffusivity and time step values (1D walk)
P = 1/(sqrt(4*pi*D*tau))*exp(-(x.^2/(4*D*tau)));

% Define the PDF given the diffusivity and time step values (2D walk)
% P = 1/(sqrt(8*pi*D*tau))*exp(-(x.^2/(8*D*tau)));

% Verify that the PDF is defined for all necessary x values
cdf = cumtrapz(x,P);
if cdf(end) < 0.99
    error_message = strcat(['WARNING. The area under the PDF is only ' num2str(cdf(end)) '.']);
    disp(error_message)
end

end

% plot(x,P)
% title('PDF')
% xlabel('x')
% ylabel('PDF')