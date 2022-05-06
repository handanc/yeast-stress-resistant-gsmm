%% Load Model Files
clear;clc;
load ale_models
experiments = { 'ethanolb2','ethanolb8', 'caffeine', 'coniferylaldehyde', 'iron', 'nickel', 'phenylethanol', 'silver', 'wildtype'};
figlegend = {'B2 (ethanol)', 'B8 (ethanol)', 'CAF905-2 (caffeine)', 'BH-13 (coniferyl aldehyde)', 'M8FE (iron)', 'M9 (nickel)', 'C9 (phenylethanol)', '2E (silver)', 'Reference'}

 
%% Run Iterative FBA for ATP minimization while a certain level of growth is required
rxn_biomass = 753;
rxn_atp = find(contains(ale_models.wildtype.rxns, 'arm_r_0226')); 

for exp = 1:length(experiments) 
    experiment = experiments{exp} ; 
    % Set model
    tmp = ale_models.(experiment);
    %tmp = changeRxnBounds(tmp, exchange_rxns_uptakes, 0, 'b');
    %tmp = changeRxnBounds(tmp, exchange_rxns_uptakes_allowed, Inf, 'u');
    tmp = changeRxnBounds(tmp, 'r_1714_REV', 10, 'b'); % glucose = 10
    
    disp(sprintf('Solving iterative FBA for %s.', experiment))
    % Get Biomass Linspace
    b = linspace(0.5,0,50);
    biopt2.(experiment) = table('Size', [length(b),2], 'VariableTypes', ["double", "double"], 'VariableNames', ["Growth", "ATP"]);
    
    % Optimize for ATP M
    tmp.c = zeros(length(tmp.rxns) ,1);
    tmp.c(rxn_atp) = 1 ; % change objective
    tmp = changeRxnBounds(tmp, tmp.rxns(rxn_biomass), 1000, 'u');  
    for i = 1:length(b)
        bm = b(i);
        tmp = changeRxnBounds(tmp, tmp.rxns(rxn_biomass), bm, 'l');  

        solatp = optimizeCbModel(tmp, 'min');
        
        if ~isempty(solatp.x) 
            biopt2.(experiment).Growth(i) = solatp.x(rxn_biomass);
            biopt2.(experiment).ATP(i) = solatp.f;
            biopt2_solutions.(experiment)(:,i) = solatp.x;
        else
            biopt2.(experiment).Growth(i) = nan;
            biopt2.(experiment).ATP(i) = nan; 
            biopt2_solutions.(experiment)(:,i)  = zeros(length(tmp.rxns) ,1);
        end
       
    end
end

%% Figures
load randomsampling_90th_5k


 
figure('Color',[1 1 1], 'Position', [0 0 2000 800]); 
    subplot(1,2,1)
for exp = 1:length(experiments) 
     ggplot();
    experiment = experiments{exp} ; 
    res = biopt2.(experiment);
    plot(res.ATP, res.Growth, 'Color', ale_models.(experiment).color, 'LineWidth', 1.5); hold on;
    xlabel('ATP Synthase (mmol/gDWh)' )
    ylabel('Growth Rate (1/h)')
end
    legend(figlegend, 'Location', 'northwest', 'box', 'on', 'FontSize', 12,'NumColumns', 2)
 
    subplot(1,2,2)
for exp = 1:length(experiments) 
    ggplot();
    experiment = experiments{exp} ;  
    points = solutions90.(experiment);
    scatter(points(rxn_atp,:), points(rxn_biomass,:), 55, ale_models.(experiment).color, 'x'); hold on;
end
    xlabel('ATP Synthase (mmol/gDWh)' )
    ylabel('Growth Rate (1/h)')
    legend(figlegend, 'Location', 'northeast', 'box', 'on', 'FontSize', 12)
    
print('figures_pdf/f10_2_biopt_atp_zoomed','-dpdf', '-r600')

% Zoomed section
% figure('Color',[1 1 1], 'Position', [0 0 900 800]);
% subplot(2,2,4)
% for exp = 1:length(experiments) 
%      ggplot();
%     experiment = experiments{exp} ; 
%     res = biopt2.(experiment);
%     plot(res.ATP, res.Growth, 'Color', ale_models.(experiment).color, 'LineWidth', 1.5); hold on;
% end