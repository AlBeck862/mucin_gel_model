function arr_no_trailing_zeros = no_trailing_zeros(arr)
% NO_TRAILING_ZEROS Removes trailing zeros from a vector.

if ~any(any(arr)) %if the vector is all zeros, return an empty matrix
    arr_no_trailing_zeros = [];
elseif arr(end)==0 && arr(end-1)==0 && arr(end-2)==0 && arr(end-3)==0 %if at least the last four elements are zeros, search for the start of the chain of zeros and eliminate it
    idxs_zeros = strfind(arr,[0,0,0,0]);
    arr_no_trailing_zeros = arr;
    arr_no_trailing_zeros(idxs_zeros(1):end) = [];
elseif arr(end)==0 && arr(end-1)==0 && arr(end-2)==0 %if the last three elements are zeros, eliminate them
    arr_no_trailing_zeros = arr(1:end-3);
elseif arr(end)==0 && arr(end-1)==0 %if the last two elements are zeros, eliminate them
    arr_no_trailing_zeros = arr(1:end-2);
elseif arr(end)==0 %if the last element is zero, eliminate it
    arr_no_trailing_zeros = arr(1:end-1);
else
    arr_no_trailing_zeros = arr;
end

end