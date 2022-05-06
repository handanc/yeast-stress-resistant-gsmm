clear;clc;
load ale_models
experiments = { 'ethanolb2','ethanolb8', 'caffeine', 'coniferylaldehyde', 'iron', 'nickel', 'phenylethanol', 'silver', 'wildtype'};
figlegend = {'B2 (ethanol)', 'B8 (ethanol)', 'CAF905-2 (caffeine)', 'BH-13 (coniferyl aldehyde)', 'M8FE (iron)', 'M9 (nickel)', 'C9 (phenylethanol)', '2E (silver)', 'Reference'}

% Collect Exchange solutions
[~, exchange_rxns_all]  = getExchangeRxns(ale_models.wildtype);

% Uptake rxns
exchange_rxns_uptakes = ale_models.wildtype.rxns(exchange_rxns_all(contains(ale_models.wildtype.rxnNames(exchange_rxns_all), '(reversible')));
exchange_rxns_uptakes_allowed = {'ammonium exchange (reversible)' 
                                'D-glucose exchange (reversible)'
                                'H+ exchange (reversible)'       
                                'iron(2+) exchange (reversible)' 
                                'oxygen exchange (reversible)'   
                                'phosphate exchange (reversible)'
                                'potassium exchange (reversible)'
                                'sodium exchange (reversible)'   
                                'sulphate exchange (reversible)' 
                                'chloride exchange (reversible)' 
                                'Cu2(+) exchange (reversible)'   
                                'Mn(2+) exchange (reversible)'   
                                'Zn(2+) exchange (reversible)'  
                                'Mg(2+) exchange (reversible)'   
                                'Ca(2+) exchange (reversible)' };
[~, ia] = intersect(ale_models.wildtype.rxnNames, exchange_rxns_uptakes_allowed, 'stable');
exchange_rxns_uptakes_allowed = ale_models.wildtype.rxns(ia);
clear ia

%% Perform MOMA
MOMAresults.maxgrowth = table();
MOMAresults.maxgrowth.rxns = (1:length(ale_models.wildtype.rxns))'; 

MOMAresults.fix_gur_10 = table();
MOMAresults.fix_gur_10.rxns = (1:length(ale_models.wildtype.rxns))'; 

% for r_id = 2 % first one is skipped
%     
    for exp = 1:length(experiments)-1 
        experiment = experiments{exp} ;

        % Models
        modelWT = ale_models.wildtype;     % WT model
        modelWT = changeRxnBounds(modelWT, exchange_rxns_uptakes, 0, 'b');
        modelWT = changeRxnBounds(modelWT, exchange_rxns_uptakes_allowed, Inf, 'u'); 
        modelWT.ub(modelWT.ub == Inf) = 1000;
        
        modelALE = ale_models.(experiment);    % ALE model  
        modelALE = changeRxnBounds(modelALE, exchange_rxns_uptakes, 0, 'b');
        modelALE = changeRxnBounds(modelALE, exchange_rxns_uptakes_allowed, Inf, 'u'); 
        modelALE.ub(modelALE.ub == Inf) = 1000;

%         if r_id == 1 % Maximize growth
%             [solALE, solWT] = MOMA(modelWT, modelALE, 'max', 1, 1); % line 133, 'one' instead of true 
%             if exp == 1
%                 MOMAresults.maxgrowth.wildtype = solWT.x;
%             end
%                 MOMAresults.maxgrowth.(experiment) = solALE.x;
%          
%         elseif r_id == 2  % Fix glucose uptake at 10
            gup = find(contains(ale_models.wildtype.rxnNames, 'D-glucose exchange (reversible)'));
            gur = 10;
            modelWT  = changeRxnBounds(modelWT, modelWT.rxns(gup), gur, 'b');
            modelALE = changeRxnBounds(modelALE, modelALE.rxns(gup), gur, 'b');
            [solALE, solWT] = MOMA(modelWT, modelALE, 'max', 1, 1); % line 133, 'one' instead of true 
                if exp == 1 % save only once
                MOMAresults.fix_gur_10.wildtype = solWT.x;
                end
            MOMAresults.fix_gur_10.(experiment) = solALE.x;
%         end
%     end 
end; clear exp experiment gur modelALE modelWT solALE solWT grate orglb orgub gro gup atp r_id

% save('data/moma_results.mat','MOMAresults')

%% Investigate results
res = MOMAresults.fix_gur_10;
%res{:,2:end} = res{:,2:end}*1e6; % convert to nmol

% Calculate std
res.std = std(table2array(res(:,3:10)),0,2); % without wildtype

% Calculate range
res.range = range(table2array(res(:,3:10)),2);

% Cumulative Distances 
tmp = table();
for exp = 1:length(experiments)-1
    experiment = experiments{exp};
    tmp.(experiment) = abs(res.wildtype - res.(experiment));
end; res.cumDist = sum((table2array(tmp)), 2);
clear exp experiment tmp  

% Remove all zero reactions
%res(all(table2array(res(:,2:10)) == 0, 2),:) = [];

% Remove rows with (No) in their names 
% res( find(contains(ale_models.wildtype.rxnNames, ' (No')),:) = []; 

% Extract proteins
draw_pos = find( contains(ale_models.wildtype.rxns, 'draw_'), 1, 'first');
res_prot = res(res.rxns > draw_pos-1, :);
res_prot{:,2:end} = res_prot{:,2:end}*1e6; % convert to nmol
res_rxn = res(res.rxns < draw_pos, :);
clear draw_pos

%% FIGURES
fig.nbox = 10;
fig.selection = 'cumDist';
fig.title = 'Cumulative Distances from WT of Protein Usage Fluxes  (Top 50 Proteins) - GUR10';
fig.data = sortrows(res_prot, fig.selection,'descend');
fig.label = '';

fig.names = ale_models.wildtype.rxnNames(table2array(fig.data(1:fig.nbox,1)));
% FOR RXN
%fig.names(contains(fig.names, ' (arm)')) = extractBefore(fig.names(contains(fig.names, ' (arm)')), ' (arm');

% FOR PROTEIN:
fig.names = extractAfter(fig.names, 10);
[~,~,ib] = intersect(fig.names, ale_models.wildtype.enzymes, 'stable');
fig.names = ale_models.wildtype.enzNames(ib);
%%
figure('Color',[1 1 1]); ggplot(); 
bar(fig.data.(fig.selection)(1:fig.nbox), 'k', 'FaceAlpha', 0.15, 'EdgeAlpha', 0, 'BarWidth', 0.5);
%fig.line = plot(fig.data.(fig.selection)(1:fig.nbox), 'LineWidth', 2);
ylabel('Cumulative Flux Distances')
%yticks(0:250:2000)
hold on; 
yyaxis right; 
scatter(1:fig.nbox, fig.data.wildtype(1:fig.nbox), 500, '*', 'k'); hold on;

for exp = 1:length(experiments)-1
    experiment = experiments{exp};
    fig.sc = scatter(1:fig.nbox, fig.data.(experiment)(1:fig.nbox), 50, ale_models.(experiment).color, 'filled');hold on;
end; clear exp experiment

xticks(1:fig.nbox); xticklabels(fig.names); 

ylabel('Enzyme Usage Fluxes (nmol/gDWh)')
%yticks(0:150:1200)
%fig.line = bar(fig.data.(fig.selection)(1:fig.nbox), 'k', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
lgnd  = legend([{'Cumulative distance'} {'Reference'} figlegend(1:end-1)], 'box', 'on') ;
%title(fig.title)

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(lgnd,'color','white');

print('figures_pdf/f9_moma_enzymes','-dpdf', '-r600')