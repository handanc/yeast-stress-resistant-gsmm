function [tickLabels] = multiLineLabels()

row1 = {'     B2','    B8','CAF905-2','   BH-13',' M8FE','   M9','         C9','   2E','Reference',''};
row2 = {'(ethanol)','(ethanol)',' (caffeine)',' (coniferyl','  (iron)','(nickel)','(phenylethanol)','(silver)','',''};
row3 = {'','','',' aldehyde)','','','','','',''};

labelArray = [row1; row2; row3]; 
tickLabels = strtrim(sprintf('%s\\newline%s\\newline%s\n', labelArray{:}));
end