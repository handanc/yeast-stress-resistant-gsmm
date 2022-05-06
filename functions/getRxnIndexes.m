function [index] = getRxnIndexes(model)

index.growth = find(strcmpi(model.rxnNames,'growth'));
index.glu    = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'));
index.oxy    = find(strcmpi(model.rxnNames,'oxygen exchange (reversible)'));
index.co2    = find(strcmpi(model.rxnNames,'carbon dioxide exchange'));
index.eth    = find(strcmpi(model.rxnNames,'ethanol exchange'));
index.ace    = find(strcmpi(model.rxnNames,'acetate exchange'));
index.gly    = find(strcmpi(model.rxnNames,'glycerol exchange'));
index.for    = find(strcmpi(model.rxnNames,'formate exchange (reversible)'));

end
