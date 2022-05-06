
%% Table with minFVA FBA maxFBA

for exp = 1:length(experiments) 
    experiment = experiments{exp} ;
    
    % Table: MinFlux - FBA Sol - MaxFlux
    flex.(experiment) = [FVA.(experiment)(:,1:2) fba_solvectors.free.(experiment) FVA.(experiment)(:,3) ];
    
    % Check if the FBA sol is within 20% of FVA
    tmp = [(flex.(experiment)(:,2) > flex.(experiment)(:,3) * .8), (flex.(experiment)(:,4) < flex.(experiment)(:,3) * 1.2 )];
    tmp = all(tmp, 2);
    
    % Collect reactions that are in the interval
    interval.(experiment) = flex.(experiment)(tmp,:);
    
end; clear exp experiment ;

 
%% Range Interval Figures of Rxns (Not Used)
for exp = 1:length(experiments) 
    experiment = experiments{exp} ;
    figure('Color',[0.95 0.95 0.95]);
    bar(interval.(experiment)(:,4), 'EdgeColor', 'none');hold on;
    bar(interval.(experiment)(:,2), 'w', 'EdgeColor', 'none');hold on;
    n = length(interval.(experiment));
    scatter(1:n, interval.(experiment)(:,3), '*');
    xticks(1:n)
    names = ale_models.wildtype.rxnNames(interval.(experiment)(:,1));
    xticklabels(names)
    xtickangle(90)
    title(sprintf("%s (n=%d)", experiment, n))
    
end; clear exp experiment ;