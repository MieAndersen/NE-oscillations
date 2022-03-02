% function Times = TriggerPoints(data,threshold,minISI)
% returns triggered points with a inter-trigger-interval of 
% at least minISI samples
% triggered points are defined as points that exceeds 'threshold'
% see also TriggerPointsEnd
% 
% 20100601 takata norio modified. empty handling.
function Times = TriggerPoints(data,threshold,minISI)

data = data(:);
tempTimes = find(data>threshold);
if isempty(tempTimes)
    Times = tempTimes;
else
    Times = tempTimes([true;diff(tempTimes)>minISI]);
end


