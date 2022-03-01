function [onset, offset] = binary_to_OnOff(binary_vector)
% BINARY_TO_ONOFF takes a logical vector and returns onsets and offsets. If
% the binary vector begins with a 1 (true), 0 is put at beginning of
% onsets. If the binary vector ends with 1 the length of binary is put at
% the end of offsets.
% 
% Input arguments:
%   binary_vector: is a logical vector
% Output arguments:
%   onset: is a column vector of indices for transitions of 0>1
%   offset:is a column vector of indices for transitions of 1>0

    onset = find(diff(binary_vector)==1)';
    offset = find(diff(binary_vector)==-1)';
    
    if onset(1) > offset(1) % in cases where the binary vector is 1 from begnning
        onset = [0; onset];
    end
    
    if onset(end) > offset(end) % in cases where binary vector is 1 at the end
        offset = [offset; length(binary_vector)];
    end
    
end
