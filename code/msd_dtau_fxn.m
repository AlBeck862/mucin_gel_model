function f = msd_dtau_fxn(x,lagT,msd)

for i=1:1:length(lagT)
    msd_fit(i) = 4*x(1)*lagT(i)^x(2);

%     msd_fit(i) = 4*exp(-abs(x(1)))*lagT(i)^exp(-abs(x(2)));
%     msd_fit(i) = 4*10^(-abs(x(1)))*lagT(i)^10^(-abs(x(2)));    
end

% Normalize the result to improve accuracy
f = sum(((msd-msd_fit)./msd).^2);

% f = sum((log(msd+1)-log(msd_fit+1)).^2);
% f = sum((log10(msd+1)-log10(msd_fit+1)).^2);

end

