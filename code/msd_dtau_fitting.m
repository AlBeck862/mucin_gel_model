function [D,a] = msd_dtau_fitting(guess_D,guess_a,lagT,msd)

% g_D =log(1/guess_D);
% g_a =log(1/guess_a);

% g_D =log10(1/guess_D);
% g_a =log10(1/guess_a);

g_D=guess_D;
g_a=guess_a;

options=optimset('TolX',10^(-15),'TolFun',10^(-15),'MaxIter',100000,'MaxFunEvals',100000,'Display','off');
[vars] = fminsearch(@msd_dtau_fxn,[g_D,g_a],options,lagT,msd);

% vars(1)
% vars(2)

% D=exp(-abs(vars(1)));
% a=exp(-abs(vars(2)));

% D=10^(-abs(vars(1)));
% a=10^(-abs(vars(2)));

D=vars(1);
a=vars(2);

end

