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

%% Run Robustness Analysis

for exp = 1:length(experiments) 
    experiment = experiments{exp} ;
    model = ale_models.(experiment);
    
    model = changeRxnBounds(model, exchange_rxns_uptakes, 0, 'b');
    model = changeRxnBounds(model, exchange_rxns_uptakes_allowed, Inf, 'u'); 
    
    [controlFlux.(experiment).glu, objFlux.(experiment).glu] = robustnessAnalysis(model, model.rxns(id.glu), 50,0);
    [controlFlux.(experiment).oxy, objFlux.(experiment).oxy] = robustnessAnalysis(model, model.rxns(id.oxy), 50,0);

end
save('data/robustness_results.mat', 'controlFlux', 'objFlux')

%% Robustness Plot : Growth vs GUR Figure
figure('Color',[1 1 1], 'WindowState', 'Maximized');

% Background Colors
    subplot(2,2,1);   
    r1 = rectangle('Position',[0 0 50 0.5]')
    r1.FaceColor = [0, 0, 0, 0.05];
    r1.EdgeColor = [0, 0, 0, 0.05];
    r1.LineWidth = 1;
    
    subplot(2,2,2);   
    r1 = rectangle('Position',[0 0 10 0.5]')
    r1.FaceColor = [0, 0, 0, 0.05];
    r1.EdgeColor = [0, 0, 0, 0.05];
    r1.LineWidth = 1;
    
for exp = 1:length(experiments) 
    experiment = experiments{exp} ;
    
    subplot(2,2,1);    ggplot();
    plot(controlFlux.(experiment).glu, objFlux.(experiment).glu, '-o', 'MarkerSize',2, 'Color', ale_models.(experiment).color, 'LineWidth', 2);
    xlabel(sprintf('Glucose Uptake Rate (mmol/gDWh)', ale_models.wildtype.rxnNames{id.glu}) )
    ylabel('Growth Rate (1/h)')
    legend(figlegend, 'Location', 'northeast', 'box', 'on', 'FontSize', 12)
    
    subplot(2,2,2);    ggplot();
    plot(controlFlux.(experiment).oxy, objFlux.(experiment).oxy, '-o', 'MarkerSize',2, 'Color',ale_models.(experiment).color, 'LineWidth', 2);
    xlabel(sprintf('Oxygen Uptake Rate (mmol/gDWh)', ale_models.wildtype.rxnNames{id.oxy}) )
    ylabel('Growth Rate (1/h)')
    legend(figlegend, 'Location', 'northeast', 'box', 'on', 'FontSize', 12)
    
    subplot(2,2,3);    ggplot();
    plot(controlFlux.(experiment).glu, objFlux.(experiment).glu, '-o', 'MarkerSize',2, 'Color', ale_models.(experiment).color, 'LineWidth', 2);
    xlabel(sprintf('Glucose Uptake Rate (mmol/gDWh)', ale_models.wildtype.rxnNames{id.glu}) )
    xlim([0 50])
    ylabel('Growth Rate (1/h)')
    
    subplot(2,2,4);    ggplot();
    plot(controlFlux.(experiment).oxy, objFlux.(experiment).oxy, '-o', 'MarkerSize',2, 'Color',ale_models.(experiment).color, 'LineWidth', 2);
    xlabel(sprintf('Oxygen Uptake Rate (mmol/gDWh)', ale_models.wildtype.rxnNames{id.oxy}) )
    xlim([0 10])
    ylabel('Growth Rate (1/h)')

hold on;
end

%print('figures_pdf/f4_robustness_analysis','-dpdf', '-r600')


%% Run FBA

fba_solvectors.free.all(:, exp+1) = (1:length(ale_models.wildtype.rxns))';
    
for exp = 1:length(experiments) 
    experiment = experiments{exp} ;
    fprintf('\nSolving FBA for the model: %s \n', experiment)
    tmp = ale_models.(experiment);
    tmp = changeRxnBounds(tmp, exchange_rxns_uptakes, 0, 'b');
    tmp = changeRxnBounds(tmp, exchange_rxns_uptakes_allowed, Inf, 'u'); 
    
    % FBA with no constraints
    sol = optimizeCbModel(tmp, 'max');
    fba_solvectors.free.(experiment) = sol.x;
    fba_solvectors.free.all(:, exp+1) = sol.x;
    fba_shadowprices.free.rxns.(experiment) = sol.w;
    fba_shadowprices.free.mets.(experiment) = sol.y;
    fba_shadowprices.free.slacks.(experiment) = sol.s;
    
    % FBA with GUR iteration
    for gur = 1:30
        tmp_name = sprintf('gur%d',gur);
        tmp = changeRxnBounds(tmp, 'r_1714_REV', gur, 'b');
        sol = optimizeCbModel(tmp, 'max');
        fba_solvectors.(tmp_name).(experiment) = (1:length(ale_models.wildtype.rxns))';
        fba_solvectors.(tmp_name).(experiment) = sol.x;
        if min(sol.x) < 0
            fprintf('OUTSIDE OF BOUNDARIES: %d ', gur)
        end
        if gur == 10
         fba_shadowprices.gur10.rxns.(experiment) = sol.w;
         fba_shadowprices.gur10.mets.(experiment) = sol.y;
         fba_shadowprices.gur10.slacks.(experiment) = sol.s;
        end
    end
end; clear exp experiment tmp sol gur name;

%save('data/fba_results.mat', 'fba_solvectors', 'fba_shadowprices') 

%% FBA for isoenzyme check

sol = table();
sol.n = (1:length(ale_models.ethanolb2.rxnNames))';
sol.rxns = ale_models.ethanolb2.rxnNames;
for exp = 1%:length(experiments) 
    experiment = experiments{exp} ;
    fprintf('\nSolving FBA for the model: %s \n', experiment)
    tmp = ale_models.(experiment);
    tmp = changeRxnBounds(tmp, exchange_rxns_uptakes, 0, 'b');
    tmp = changeRxnBounds(tmp, exchange_rxns_uptakes_allowed, Inf, 'u'); 
    
    % FBA 
    tmp = changeRxnBounds(tmp, 'r_1714_REV', 10, 'b');
    tmp = changeRxnBounds(tmp, tmp.rxns(7184), 0, 'b');
    tmpsol = optimizeCbModel(tmp, 'max');
    tmpsol.f
    sol.fba = tmpsol.x;
   
end

tdh = ale_models.ethanolb2.enzymes(contains(ale_models.ethanolb2.enzNames, 'TDH'));
sol(find(contains(sol.rxns, tdh)),:)



 

%% Find Deactivated and Activated Rxns
for exp = 1:length(experiments)-1
    experiment = experiments{exp} ;
    
    % For free uptake
        tmp_solwt = fba_solvectors.free.wildtype;
        tmp_solale = fba_solvectors.free.(experiment);
        
        % Find passive and active rxns in wt
        tmp_zeros    = find(tmp_solwt == 0); % n = 7225
        tmp_nonzeros = find(tmp_solwt > 0); % n = 912
        
        % Compare them to evolved model
        iszero = [tmp_zeros tmp_solale(tmp_zeros)];
        iszero(iszero(:,2) == 0,:) = []; % remove nonzeros bc no change there
        
        isnonzero = [tmp_nonzeros tmp_solale(tmp_nonzeros)];
        isnonzero(isnonzero(:,2) > 0,:) = []; % remove zeros bc no change there
        
        activated_rxns = iszero(:,1); % n = 66
        deactivated_rxns = isnonzero(:,1); % n = 55
        
        % Store the results as: RxnIDS - WT Flux - EV Flux
        rxns_activated.free.(experiment)   = [activated_rxns tmp_solwt(activated_rxns) tmp_solale(activated_rxns)];
        rxns_deactivated.free.(experiment) = [deactivated_rxns tmp_solwt(deactivated_rxns) tmp_solale(deactivated_rxns)];
 
    % For GUR10 
        tmp_solwt = fba_solvectors.gur10.wildtype;
        tmp_solale = fba_solvectors.gur10.(experiment);
        
        % Find passive and active rxns in wt
        tmp_zeros    = find(tmp_solwt == 0); % n = 7225
        tmp_nonzeros = find(tmp_solwt > 0); % n = 912
        
        % Compare them to evolved model
        iszero = [tmp_zeros tmp_solale(tmp_zeros)];
        iszero(iszero(:,2) == 0,:) = []; % remove nonzeros bc no change there
        
        isnonzero = [tmp_nonzeros tmp_solale(tmp_nonzeros)];
        isnonzero(isnonzero(:,2) > 0,:) = []; % remove zeros bc no change there
        
        activated_rxns = iszero(:,1); % n = 66
        deactivated_rxns = isnonzero(:,1); % n = 55
        
        % Store the results as: RxnIDS - WT Flux - EV Flux
        rxns_activated.gur10.(experiment)   = [activated_rxns tmp_solwt(activated_rxns) tmp_solale(activated_rxns)];
        rxns_deactivated.gur10.(experiment) = [deactivated_rxns tmp_solwt(deactivated_rxns) tmp_solale(deactivated_rxns)];

end;
%save('data/fba_activated_deactivated_rxns.mat', 'rxns_deactivated', 'rxns_activated')


%% ATP - Growth - GUR Figure

% Collect fluxes in a single vector
model = ale_models.wildtype;

rxn_biomass = 753;
metATP = {'s_0434', 's_0435', 's_0437', 's_0438', 's_0439'};
rxn_atp = findRxnsFromMets(model, metATP, false, 'producersOnly', true);
rxn_atp = find(contains(ale_models.wildtype.rxns, rxn_atp)); 
rxn_glu = find(contains(ale_models.wildtype.rxns, 'r_1714_REV')); 

fluxes = table();
for exp = 1:length(experiments)
    experiment = experiments{exp} ; 
    fluxes.name(exp) = {ale_models.(experiment).legendName};
    fluxes.gro(exp) = fba_solvectors.free.(experiment)(rxn_biomass);
    fluxes.glu(exp) = fba_solvectors.free.(experiment)(rxn_glu);
    fluxes.atp(exp) = sum(fba_solvectors.free.(experiment)(rxn_atp));
end

%% Figure
fluxes = sortrows(fluxes,'gro','descend');
figure('Color',[1 1 1], 'WindowState', 'maximized');
ggplot();

x = 1:length(experiments);

bar(fluxes.gro, 'FaceAlpha', 0.3, 'FaceColor', 'k', 'EdgeColor', 'none', 'BarWidth', 0.4)
ylabel('Growth Rate (1/h)')
xticks(1:length(experiments))
label = extractBefore(fluxes.name, ' ');
label(find(cellfun(@isempty,label))) = {'Reference'};
xticklabels(label)
xtickangle(0)

yyaxis right
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

scatter(x,fluxes.glu,50,'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.85,0.33,0.10])
hold on;
scatter(x,fluxes.atp,50,'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.00,0.45,0.74])
lsline()
ylabel('Fluxes (mmol/gDWh)')

legend({'Growth Rates', 'GUR','ATP', },'AutoUpdate','off')

print('figures_pdf/f6_fba_gro_gur_atp','-dpdf', '-r600')

%% Shadow Prices

model = ale_models.wildtype;

TableRes = table();
TableRes.rxnID = (1:length(model.rxns))';
TableRes.rxnNames = model.rxnNames;

for exp = 1:length(experiments) 
    experiment = experiments{exp} ; 
    TableRes.(experiment) = fba_shadowprices.gur10.rxns.(experiment);
end
TableRes(all(TableRes{:,3:11} == 0, 2),:) = [];
TableRes.std = std(TableRes{:,3:11},0, 2);

%TableRes(contains(TableRes.rxnNames, 'exchange'),:) = [];
%TableRes(contains(TableRes.rxnNames, '(No'),:) = [];
TableRes(any(TableRes{:, 3:11} > -1, 2),:) = []; % -1 is the growth, remove less effectives
TableRes = sortrows(TableRes,'ethanolb2','ascend');

revid = contains(TableRes.rxnNames, ' (rev'); % not required in figure
TableRes.rxnNames(revid) = extractBefore(TableRes.rxnNames(revid), ' (rev');

exvid = contains(TableRes.rxnNames, ' exchange'); % not required in figure
TableRes.rxnNames(exvid) = extractBefore(TableRes.rxnNames(exvid), ' exchange');

% HARD CODED <<<
TableRes.rxnNames{6} = 'protein pool limitation' ;
TableRes.rxnNames{10} = 'biomass production';

TableRes.mean = mean(TableRes{:,3:10},2);

%  Results Table
for exp = 1:length(experiments)
experiment = experiments{exp} ; 
TableRes.(experiment) = (TableRes.(experiment) ./ TableRes.wildtype)
end

%% Figure
row1 = {'biotin','riboflavin','FMN','5-formyltetra','4-amino',' protein ','ergosterol','       ergosta-','zymosterol','  biomass'};
row2 = {'','','','  hydrofolic','benzoate','   pool','','  5,7,22,24(28)','','production'};
row3 = {'','','','      acid','','limitation','','tetraen-3-beta-ol','',''};
labelArray = [row1; row2 ;row3]; 
tickLabels = strtrim(sprintf('%s\\newline%s\\newline%s\n', labelArray{:}));

figure('Color',[1 1 1]); 
subplot(3,1,1:2);ggplot();
   H=bar(TableRes{:,3:10})
   xticks(1:length(TableRes.rxnNames))
   xticklabels(tickLabels)
  % xtickangle(90)
   ylabel('Shadow Prices (ale/wt)')
for exp = 1:length(experiments)-1
    experiment = experiments{exp} ; 
       H(exp).FaceColor = 'flat';
       H(exp).CData = ale_models.(experiment).color;
end
  legend(figlegend, 'Location','northoutside', 'NumColumns', 4, 'box', 'on')
    
print('figures_pdf/f11_fba_shadowprices','-dpdf', '-r600')
