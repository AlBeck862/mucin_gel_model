function arr_no_trailing_zeros = no_trailing_zeros(arr)
% NO_TRAILING_ZEROS Removes trailing zeros from a vector.

if ~any(any(arr)) %if the vector is all zeros, return an empty matrix
    arr_no_trailing_zeros = [];
else %only keep values prior to the trailing zeros
    last_nonzero_value_idx = find(arr,1,'last');
    arr_no_trailing_zeros = arr(1:last_nonzero_value_idx);
end

end