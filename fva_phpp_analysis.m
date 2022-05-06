%% Load Models
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

id = getRxnIndexes(ale_models.wildtype);

%% Run FVA on ALE
for exp = 1:length(experiments)
    experiment = experiments{exp};

    tmp = ale_models.(experiment);
    tmp = changeRxnBounds(tmp, exchange_rxns_uptakes, 0, 'b');
    tmp = changeRxnBounds(tmp, exchange_rxns_uptakes_allowed, Inf, 'u'); 
    tmp.ub(tmp.ub == Inf) = 1000;

    disp(sprintf("Solving for: %s", experiment));
    [minFlux, maxFlux]  = fluxVariability(tmp,  90, 'max', tmp.rxns);
    FVA.(experiment) = [(1:length(minFlux))' minFlux maxFlux];
    clear minFlux maxFlux
end; clear tmp experiment


%% Analyze the results

% Precreate tables
enzyme_ranges = table();
enzyme_ranges.RxnIds = FVA.wildtype(:,1);
rxn_ranges = table();
rxn_ranges.RxnIds = FVA.wildtype(:,1);

% Calculate ranges
for exp = 1:length(experiments)
    experiment = experiments{exp};
    fvares = FVA.(experiment);
    fvares(fvares(:,3) > 1000,3) = 1000; % Convert Infs to 1000 for calculation
    rxn_ranges.(experiment)    =  fvares(:,3) - fvares(:,2);
    enzyme_ranges.(experiment) = (fvares(:,3) - fvares(:,2))*1e+6; % converted to nmol bc numbers are too small for enzymes
end; clear exp fvares

% Divide results for enzyme and metabolic rxns
draw_pos = find( contains(ale_models.wildtype.rxns, 'draw_'), 1, 'first'); % position of the first enzyme rxn

enzyme_ranges(enzyme_ranges{:,1} < draw_pos,:) = []; % Obtain only proteins
enzyme_ranges(end,:) = []; % Remove the last rxn: prot_pool

rxn_ranges(rxn_ranges{:,1} >= draw_pos,:) = []; % Obtain only rxns

% Calculate standart deviations
enzyme_ranges.std = std(enzyme_ranges{:,2:length(experiments)+1},0,2); 
rxn_ranges.std = std(rxn_ranges{:,2:length(experiments)+1},0,2); 


% Insert Protein Names
prot_Names = extractAfter(ale_models.wildtype.rxnNames(enzyme_ranges.RxnIds), 10);
[~, ia, ib] = intersect(prot_Names, ale_models.wildtype.enzymes, 'stable');
prot_Enzymes = ale_models.wildtype.enzNames(ib);
prot_Genes = ale_models.wildtype.enzGenes(ib);
enzyme_ranges = addvars(enzyme_ranges, prot_Names, prot_Enzymes, prot_Genes, 'After', 'RxnIds');

rxn_Names = ale_models.wildtype.rxnNames(rxn_ranges.RxnIds);
rxn_ranges = addvars(rxn_ranges, rxn_Names, 'After', 'RxnIds');

% Sort by Standart Devs.
enzyme_ranges = sortrows(enzyme_ranges, 'std', 'descend');
rxn_ranges = sortrows(rxn_ranges,  'std', 'descend'); 

% Save data
% save('data/fva_results.mat','FVA', 'enzyme_ranges', 'rxn_ranges')


%% Obtain the table for the most divergent enzymes by std.
fva_ranked_enzymes = enzyme_ranges(1:50,:); % get the first 50 
fva_ranked_enzymes.range = max(fva_ranked_enzymes{:,5:13}, [],2) - min(fva_ranked_enzymes{:,5:13}, [],2);
fva_ranked_enzymes = sortrows(fva_ranked_enzymes, 'range', 'descend');
fva_ranked_enzymes = sortrows(fva_ranked_enzymes, 'std', 'descend'); 
 

%% Figure of Enzyme Ranges
enzyme_ranges = sortrows(enzyme_ranges, 'wildtype', 'descend');

figure('Color',[1 1 1], 'WindowState', 'Maximized'); 
subplot(5,1,1:3); ggplot();
    for exp = 1:length(experiments)
        experiment = experiments{exp};
        plot(enzyme_ranges.(experiment), 'LineWidth', 2, 'Color', ale_models.(experiment).color);hold on
    end
     grid on
    legend(figlegend, 'box', 'off', 'FontSize', 16)
    prot_names = extractAfter(ale_models.wildtype.rxnNames(enzyme_ranges{:,1}), 10);
    [~,~,ib] = intersect(prot_names, ale_models.wildtype.enzymes, 'stable');
    prot_names = ale_models.wildtype.enzNames(ib);
    xlim([0 length(prot_names)])
    xticks(0:50:length(prot_names))
    ylabel('Flux Ranges (nmol/gDWh)')
subplot(5,1,4:5);ggplot();
        bar(enzyme_ranges.std, 2, 'k', 'FaceAlpha', 0.8);
            grid off
    xlim([0 length(prot_names)])
    xticks(0:50:length(prot_names))
    xlabel('Enzyme Reaction #')
    ylabel('Standart Deviations')
    
    
%print('figures_pdf/fva_analysis_enzymes','-dpdf', '-r600')

%% Phenotype Phase Plane

for exp = 1:length(experiments) 
    experiment = experiments{exp} ; 
    model = ale_models.(experiment) ;
    [PPPgRates.(experiment), PPPsPrices1.(experiment), PPPsPrices2.(experiment)] = phenotypePhasePlane(model, model.rxns(id.glu), model.rxns(id.oxy), 50, 50, 50);
    %[controlFlux1, controlFlux2, objFlux] = doubleRobustnessAnalysis(model, model.rxns(ind.glu), model.rxns(ind.oxy), 20, 1, model.rxns(find(model.c)), 'max');
end
 
% save('data/ppp_results.mat','PPPgRates', 'PPPsPrices1', 'PPPsPrices2')
    

%% PhPP Figures: Countour plots
    n = 50;
    si = 1;
    figure('Color',[1 1 1])%, 'Position',  [0, 0, 400, 1080]);
for exp = 1:length(experiments)
    experiment = experiments{exp} ; 
    model = ale_models.(experiment) ;
    subplot(3,3,si); 
    ggplot;
    grid off;
    p = contour(1:n, 1:n, PPPgRates.(experiment), 'ShowText','on',  'LabelSpacing', 400);
    %view([-190 120 120])
    %p(1).EdgeColor = 'none';
    xticks(0:10:50)
    yticks(0:10:50)
    xlim([0 50])
    xlabel('Glucose Uptake Rate (mmol/gDWh)') %strrep(strcat(model.rxnNames(ind.glu),' (mmol/g DW-hr)'),'_','\_'));
    ylabel('Oxygen Uptake Rate (mmol/gDWh)') %strrep(strcat(model.rxnNames(ind.oxy),' (mmol/g DW-hr)'),'_','\_'));
    zlabel('Growth Rate (1/hr)');
    hold on;
    [row, col] = find(ismember(PPPgRates.(experiment), max(PPPgRates.(experiment)(:))));
    plot3(col,row,PPPgRates.(experiment)(row,col),'rx'); 
    hold on;
   text(col-2,row+2,string(round(max(max(PPPgRates.(experiment))),2)), 'Color', 'r');
    title(model.legendName)
    si = si +1;
end
print('figures_pdf/f5_contourplots','-dpdf', '-r600')

%% GUR O2 Ranges Table
gurTable = table();
gurTable.Names = figlegend';
gurTable.GurRange = experiments';
gurTable.GurDiff = experiments';
gurTable.O2Range = experiments';
gurTable.O2Diff = experiments';

for exp = 1:length(experiments)
    experiment = experiments{exp} ; 
    data = PPPgRates.(experiment);
    data ( data < 0.3) = 0; 
    
    gurTable.MaxGrowth(exp) = max(max(data));
    gurTable.GurRange(exp) = {sprintf('%d - %d', find(sum(data) ~= 0, 1, 'first'), find(sum(data) ~= 0, 1, 'last'))};
    gurTable.GurDiff(exp) = {find(sum(data) ~= 0, 1, 'last') - find(sum(data) ~= 0, 1, 'first') };
    gurTable.O2Range(exp) =  {sprintf('%d - %d', find(sum(data,2) ~= 0, 1, 'first'), find(sum(data,2) ~= 0, 1, 'last'))};
    gurTable.O2Diff(exp) = {find(sum(data,2) ~= 0, 1, 'last') - find(sum(data,2) ~= 0, 1, 'first') };
end
 

