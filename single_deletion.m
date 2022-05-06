%% Load Model Files
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
%% ALE
% IMPORTANT: Line 44 of singleRxnDeletion is edited to remove 'one'
% parameter in the optimization.
for exp = 1:length(experiments)
    experiment = experiments{exp} ;
    model = ale_models.(experiment);
    model = changeRxnBounds(model, exchange_rxns_uptakes, 0, 'b');
    model = changeRxnBounds(model, exchange_rxns_uptakes_allowed, Inf, 'u'); 
    [sdresults.(experiment).grRatio, sdresults.(experiment).grRateKO, sdresults.(experiment).grRateWT, sdresults.(experiment).hasEffect, sdresults.(experiment).delRxn, sdresults.(experiment).fluxSolution] = singleRxnDeletion(model, 'FBA')
end


%% Collect Results
for exp = 1:length(experiments)
    experiment = experiments{exp} ;
    colresults.grRatio(:,exp) = sdresults.(experiment).grRatio;
    colresults.grRateKO(:,exp) = sdresults.(experiment).grRateKO;
    colresults.grRateWT(:,exp) = sdresults.(experiment).grRateWT;
    colresults.hasEffect(:,exp) = sdresults.(experiment).hasEffect;
    colresults.delRxn(:,exp) = sdresults.(experiment).delRxn;

end
 colresults.experiments = experiments;
  
save('data/singledeletion_results.mat', 'sdresults', 'colresults', '-v7.3')
    
%% Analyze Data

draw_pos = find( contains(ale_models.wildtype.rxns, 'draw_'), 1, 'first');

% Metabolic Reactions
srv_rxn = table();
srv_rxn.rxns = findRxnIDs(ale_models.wildtype, colresults.delRxn(:,1));
srv_rxn.grRatio = colresults.grRatio;

srv_rxn(any(isnan(srv_rxn.grRatio), 2),:) = []; % remove unfeasible solutions if any
srv_rxn((srv_rxn.rxns > draw_pos-1),:) = []; % remove proteins
srv_rxn(all(srv_rxn.grRatio == srv_rxn.grRatio(:,1), 2), :) = []; %% remove if all the same

srv_rxn.min = min(srv_rxn.grRatio,[],2);
srv_rxn = sortrows(srv_rxn, 'min','ascend'); 
srv_rxn((srv_rxn.min > 0.95),:) = []; % remove less effectives

srv_rxn(contains(ale_models.wildtype.rxnNames(srv_rxn.rxns), ' (No'),:) = []; %% remove isozymes


% Enzymes

srv_prot = table();
srv_prot.rxns = findRxnIDs(ale_models.wildtype, colresults.delRxn(:,1));
srv_prot.grRatio = colresults.grRatio;

srv_prot((srv_prot.rxns < draw_pos),:) = []; % remove metabolic rxns
srv_prot(any(isnan(srv_prot.grRatio), 2),:) = [];% remove unfeasible solutions if any
srv_prot(all(srv_prot.grRatio == srv_prot.grRatio(:,1), 2), :) = []; %% remove all the same

srv_prot.min = min(srv_prot.grRatio,[],2);
srv_prot = sortrows(srv_prot, 'min','ascend'); 
srv_prot((srv_prot.min > 0.95),:) = []; %% remove isozymes
 
 
%% Figure
figure('Color',[1 1 1]);

    n = length(srv_prot.min);
    h = heatmap(srv_prot.grRatio(1:n,:), 'CellLabelColor', 'none',...
       'FontName', 'Helvetica', 'FontSize', 13)
   
    prot_names = extractAfter(ale_models.wildtype.rxnNames(srv_prot.rxns(1:n)), 10);
    [~,~,ib] = intersect(prot_names, ale_models.wildtype.enzymes, 'stable');
    prot_names = ale_models.wildtype.enzNames(ib);
    
    %prot_names = predef_prot_names; % names are optained from yeastgenome
    h.XDisplayLabels =(figlegend);
    h.YDisplayLabels = (prot_names) ;
    h.ColorLimits = [0 0.95];
    s = struct(h);
    s.XAxis.TickLabelRotation = 90;  % vertical

print('figures_pdf/f12_singledeletion_proteins','-dpdf', '-r600')
