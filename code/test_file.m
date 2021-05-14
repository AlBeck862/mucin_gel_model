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

[x,pdf,cdf] = gen_PDF(100,0.1);
plot(x,pdf)
title('PDF')
xlabel('Change to the step to be taken')
ylabel('Probability')

figure()
plot(x,cdf)

cdf = round(cdf,5);
% figure()
% plot(x,cdf)

displmnts = zeros(1000000,1);
for i = 1:1000000
    val = rand();
    displmnt_idx = find(cdf==round(val,5),1);
    displmnt = x(displmnt_idx);
    try
        displmnts(i)=displmnt;
    catch
        displmnts(i)=0;
    end
end

displmnts
displmnts(displmnts==0) = [];

figure()
histogram(displmnts,100)

