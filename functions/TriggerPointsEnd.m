% function Times = TriggerPointsEnd(data,threshold,minISI)
% returns triggered points with a inter-trigger-interval of 
% at least minISI samples
% triggered points are defined as "end" points that exceeds 'threshold'
% as opposed to TriggerPoints, which detects the onset of the trigger
% see also TriggerPoints
% 
% 20100601 takata norio modified. empty handling.
function Times = TriggerPointsEnd(data,threshold,minISI)

data = data(:);
tempTimes = find(data>threshold);
if isempty(tempTimes)
    Times = tempTimes;
else
    Times = tempTimes([diff(tempTimes)>minISI;true]);
end

