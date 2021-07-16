function dispVar = get_dispmnt_variation(diffusivity,x,all_cdf,diffusivities)
% GET_DISPMNT_VARIATION Return the displacement to be undergone
% by a simulated particle given a diffusivity (D) and a time lag (tau).

% Fetch the appropriate CDF for the given diffusivity

% tic
persistent cdf previous_diff

if isempty(cdf) || diffusivity ~= previous_diff
    cdf = all_cdf(diffusivities==diffusivity,:);
    cdf = round(cdf,5);
    previous_diff = diffusivity;
end
% toc

is_valid_x = 0; %used to fetch a new x value if none was found in the rounded CDF
while is_valid_x == 0
    val = rand(); %generate a random number in the range 0 to 1
    displmnt_idx = find(cdf==round(val,5),1); %get index of CDF corresponding to that random number
    displmnt = x(displmnt_idx); %get the extra displacement corresponding to that index
    
    % If no extra displacement value was found in the CDF (due to rounding), the while loop will be repeated.
    if ~isempty(displmnt)
        is_valid_x = 1;
    end
end

% Define the return variable
dispVar = displmnt;

end