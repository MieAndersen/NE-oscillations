% function EventPair = PairEvents(Beginnings,Ends)
% Given event pairs Beginnings and Ends (both a single column vectors),
% the function pairs the events with exception handling.
% Beginnings and Ends has to be in ascending order and they should have the
% same length or their length may differ by 1. (where the first event
% detection starts with an end or the last event detection does not have an
% end; cf. square wave event detection.)
%
% e.g. Beginnings = [1, 3, 6, 10]';
%      Ends = [1.1,3.5,7, 12.5]';
%      EventPair = [[1,1.1];[3,3.5];[6,7];[10,12.5]]
% e.g. Beginnings = [3, 6, 10]';
%      Ends = [1.1,3.5,7, 12.5]';
%      EventPair = [[3,3.5];[6,7];[10,12.5]]
% e.g. Beginnings = [1, 3, 6, 10]';
%      Ends = [1.1,3.5,7]';
%      EventPair = [[1,1.1];[3,3.5];[6,7]]
%
% error-prone code, but it works if event detection is performed with a
% threshold-based algorithm on a continuous signal.
% [as above, it is based on many assumptions (time order, alternating
% begin/end timings, etc]. Error-proof code to be written later, if necessary]
% 2019/03/27 Haj Hirase

function EventPair = PairEvents(Beginnings,Ends)

% make input vectors column vectors
Beginnings = Beginnings(:);
Ends = Ends(:);

if abs(length(Beginnings) - length(Ends)) > 1
    error('in appropriate vector length')
end


if length(Beginnings) == length(Ends)
    if Beginnings(1) < Ends(1)
        % case1 begin end ... begin end
        EventPair = [Beginnings(:),Ends(:)];
    else
        % case2 end begin ... end begin
        EventPair = [Beginnings(2:end),Ends(1:end-1)];
    end
else
    if length(Beginnings) > length(Ends)
        % case 3 begin end ... begin
        EventPair = [Beginnings(1:end-1),Ends(:)];
    else
        % case 4 end begin ... end
        EventPair = [Beginnings(1:end),Ends(2:end)];
    end
end