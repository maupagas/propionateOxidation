function [tickLabel] = rationalTickLabels(tickvals, denominator, integerFlag)

if integerFlag == 1,
    idInteger = mod(tickvals,1);
else
    idInteger = ones(size(tickvals));
end

for i=1:length(tickvals)
    if idInteger(i) == 0,
        tickLabel{i} = sprintf('%d',(tickvals(i)));
    else
        tickLabel{i} = sprintf('%.f/%d',tickvals(i)*denominator, denominator);
    end
end