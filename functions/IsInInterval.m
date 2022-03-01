function [OutputBool, OutputIdx, OutputData] = IsInInterval(Data,Interval)
% function [OutputBool, OutputIdx, OutputData] = IsInInterval(Data,Interval)
% It checks the if the one-dimensional  Data is in Nx2 Interval(s)
% [begin,end] for intervals defined as [begin, end) (i.e. begin
% inclusive and end exclusive: begin <= x < end)
% OutputBool is a boolean (true, false) array corresponding to each element
% of data and indicates if the element belongs to any of the defined
% intervals.
% OutputIdx is the indices to Data elements that belong to the defined
% interval(s).
% OutputData is the Data elements that belong to the defined interval(s).
% 20190402 Hajime Hirase


OutputBool = logical(zeros(length(Data),1)); 
for ii=1:length(Data)
    DataElement = Data(ii);
    for jj=1:size(Interval,1)
        if DataElement>=Interval(jj,1) && DataElement <Interval(jj,2)
            OutputBool(ii)=true;
        end
    end
end

OutputIdx = find(OutputBool);
OutputData = Data(OutputIdx);
