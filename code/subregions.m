function num_regions = subregions(het, area)
% SUBREGIONS Randomly generate a number of lattice subregions given a
% heterogeneity parameter.

max_subregions = round(1e-2*sqrt(area)); %this might need to be tuned or computed differently altogether

num_regions = 1;
for i = 1:max_subregions-1
    if rand() <= het
        num_regions = num_regions + 1;
    end
end

end