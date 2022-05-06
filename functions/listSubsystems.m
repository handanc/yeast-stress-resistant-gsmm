function [subSystems, struct] = listSubsystems(model, namelist)

  %%
draw_pos = find( contains(model.rxns, 'draw_'), 1, 'first');
subSystems = model.subSystems(1:draw_pos-1);
for i = 1:length(subSystems)
    if length(subSystems{i}) > 1
        str = '';
        for j=1:length(subSystems{i})
            str = strcat(str, {', '}, subSystems{i,1}{1,j});
        end
        subSystems(i) = str;
    else
        subSystems{i} = subSystems{i,1}{1,1};
    end
end

 %%
for i= 1:length(namelist)
    sub = namelist{i};
    struct.(sub) = find(contains(subSystems, sub));
end

end