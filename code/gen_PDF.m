function [x,pdf,cdf] = gen_PDF(D,tau)
% GEN_PDF Generate a PDF and CDF for a given diffusivity (D) and lag time
% (tau).

% Define a sufficiently broad range of x values
x = linspace(-125,125,1250000);

% Define the PDF given the diffusivity and time step values
P = 1/(sqrt(4*pi*D*tau))*exp(-(x.^2/(4*D*tau)));

% plot(x,P)
% title('PDF')
% xlabel('x')
% ylabel('PDF')

% Verify that the PDF is defined for all necessary x values
cdf = cumtrapz(x,P);
if cdf(end) < 0.975
    error_message = strcat(['ERROR. The area under the PDF is only ' num2str(cdf(end)) '.']);
    disp(error_message)
else
    pdf = P;
end

end