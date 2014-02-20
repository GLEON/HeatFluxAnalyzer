function [ meanout ] = nanmean( inputdata )

nanI = isnan(inputdata);
if sum(nanI)<length(inputdata)
    meanout = mean(inputdata(~nanI));
else
    meanout = NaN;
end


end

