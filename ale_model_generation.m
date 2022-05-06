clear;clc;
changeCobraSolver('gurobi', 'all')

%% Convert to ec model using GECKO formulation
model = readCbModel('yeastGEM.mat');
model = buildRxnGeneMat(model); % for the COBRA and RAVEN Toolbox compability
yeastGEM = model;

cd('../GECKO-2.0.0/geckomat')
tic
[ecModel,ecModel_batch] = enhanceGEM(model,'COBRA','ec_YeastGEM',1)
toc  % 40 mins

cd(currentDir)

%% Initiate Model
% load ecmodel_batch.mat
model = ecModel_batch;

experiments = { 'ethanolb2','ethanolb8', 'caffeine', 'coniferylaldehyde', 'iron', 'nickel', 'phenylethanol', 'silver'};
experimentsLegend = {'B2 (ethanol)', 'B8 (ethanol)', 'CAF905-2 (caffeine)', 'BH-13 (coniferyl aldehyde)', 'M8FE (iron)', 'M9 (nickel)', 'C9 (phenylethanol)', '2E (silver)', 'Reference'};
colors = { '#4e46f6','#1a8dec', '#37ca4e', '#fcba03', '#fc4948', '#fe7e0e', '#f429c1', '#8240e3'};
colors = hex2rgb(colors);

%% Load expressions toptables for all experiments
for i=1:length(experiments)
    toptable{i}= readExpressions(experiments{i});
end; clear i

toptable(2,:)= experiments; % Save the names to prevent confusions

% Convert Gene IDs to Systematic Names
platform = struct2cell(tdfread('differential_toptables/platform_for_genes.tsv', '\t'));
platform = cellstr(platform{1});

for exp = 1:length(experiments)
    % obtain the table for interested experiment for the rest of the script
    experiment = experiments{exp} ;
    expr.(experiment) = toptable{1, find(strcmp(toptable(2,:), experiment), 1)};
    [~, ia, ib] = intersect(expr.(experiment).Gene, platform(:,1),  'stable');
    expr.(experiment).Gene(ia) = platform(ib,2);
    expr.(experiment)(expr.(experiment).PValue > 0.05, :) = []; % second check
end; clear ia ib experiment
toptable = expr;


%% Arrange toptable experiments
fc_percentage = 100; % applied percentage
for exp = 1:length(experiments)
       
    % obtain the table for interested experiment for the rest of the script
    experiment = experiments{exp} ;
    tt_exp = toptable.(experiment);

    % Find common enzyme names
    ale_models.(experiment) = model;
    ale_models.(experiment).ale = sprintf('%s',experiment);
    ale_models.(experiment).color = colors(exp,:);

    ale_models.(experiment).enzFoldChange = zeros(length(ale_models.(experiment).enzGenes), 1);
    [~,~,ib] = intersect(ale_models.(experiment).enzGenes, table2cell(tt_exp(:,1)), 'stable');
    enzymes_toptable = table2cell(tt_exp(ib,1:2));
    ale_models.(experiment).enzFoldChangeN =  length(enzymes_toptable);
    fprintf('%d enzymes from toptable are found in %s model.\n', length(enzymes_toptable), experiment);

    for i=1:length(ale_models.(experiment).enzFoldChange)
        % find the index
        tmp = find(strcmp(ale_models.(experiment).enzGenes(i), enzymes_toptable(:,1)));
        if ~isempty(tmp) % if enzyme exists in the tt, get the fold change value (log10)
            ale_models.(experiment).enzFoldChange(i, 1) =  2^(cell2mat(enzymes_toptable(tmp, 2)) *(fc_percentage/100)) ;
        else  % if enzyme does not exist in the tt, set fold change to 1
            ale_models.(experiment).enzFoldChange(i, 1) = 1;
        end
    end

end; clear tt_exp experiment ib tmp i exp

clear expr fc_percentage enzymes_toptable colors % no more required

%% Common Genes

for exp = 1:length(experiments)
    experiment = experiments{exp} ;
    tmp = toptable.(experiment);
    commons.upregulateds.(experiment) = table2array(tmp(tmp.logFC > 1,1));
    commons.downregulateds.(experiment) = table2array(tmp(tmp.logFC < 1,1));
end
 
% Common upregulateds
tmpu = commons.upregulateds;
commons.upregulateds.list = intersect(...
                                intersect(...
                                    intersect(...
                                        intersect(...
                                            intersect(...
                                                intersect(...
                                                    intersect(tmpu.ethanolb2, tmpu.ethanolb8),...
                                                tmpu.caffeine),...
                                            tmpu.coniferylaldehyde),...
                                        tmpu.iron),...
                                    tmpu.nickel),...
                                tmpu.phenylethanol),...
                            tmpu.silver);

% Common downregulateds
tmpu = commons.downregulateds;
commons.downregulateds.list = intersect(...
                                intersect(...
                                    intersect(...
                                        intersect(...
                                            intersect(...
                                                intersect(...
                                                    intersect(tmpu.ethanolb2, tmpu.ethanolb8),...
                                                tmpu.caffeine),...
                                            tmpu.coniferylaldehyde),...
                                        tmpu.iron),...
                                    tmpu.nickel),...
                                tmpu.phenylethanol),...
                            tmpu.silver);

clear tmpu tmp

%% Apply fold changes 
for exp = 1:length(experiments)
    experiment = experiments{exp} ;
    % Change the coefficients of proteins in reactions
    fprintf('Updating coefficients for the model %s.\n', experiment);
    for i=1:length(ale_models.(experiment).enzFoldChange)
        % Find the metabolite ID of the enzyme
        enzymeName = ale_models.(experiment).enzymes{i};
        metName = ['prot_' enzymeName];
        protid = findMetIDs(ale_models.(experiment), metName);
        
        % Find the reactions that enzyme is used
        rxnsofprot = findRxnsFromMets(ale_models.(experiment), metName);
        
        % Remove draw rxn
        rxnsofprot(contains(rxnsofprot, 'draw_prot')) = [];
        rxnids = findRxnIDs(ale_models.(experiment), rxnsofprot);
        rxnName = ale_models.(experiment).rxnNames(rxnids);
        
        % Check the coefficients
        coef_prev = ale_models.(experiment).S(protid,rxnids);
        coef_fc = ale_models.(experiment).enzFoldChange(i);
        coef_new = coef_prev * coef_fc;
        
        % Update the coefficients
        ale_models.(experiment).S(protid,rxnids) = coef_new;
        
        % Save the changes to model
         ale_models.(experiment).enzCoefChanges(i,1) = {enzymeName};
         ale_models.(experiment).enzCoefChanges(i,2) = {rxnName};
         ale_models.(experiment).enzCoefChanges(i,3) = {coef_prev};
         ale_models.(experiment).enzCoefChanges(i,4) = {coef_fc};
         ale_models.(experiment).enzCoefChanges(i,5) = {coef_new};
    end
end

%% Save Models

% Add proper names for figure legends
for exp = 1:length(experiments)
    experiment = experiments{exp} ;
    ale_models.(experiment).legendName = experimentsLegend{exp};
end

% Include the wildtype model to call all of them
ale_models.wildtype = model;
ale_models.wildtype.color = hex2rgb('#88929c');
ale_models.wildtype.legendName = 'Reference';

% save('models/ale_models.mat', 'ale_models')  % done once

%%  Check expression vs. flux fold change if they are parallel  
for exp = 1:length(experiments)
    experiment = experiments{exp} ;
    if exp == 1 % only required once
        model_wt = ale_models.wildtype;
        sol_wt = optimizeCbModel(model_wt, 'max', 0, 0)
    end
    model_ale = ale_models.(experiment);
    sol_ale = optimizeCbModel(model_ale, 'max', 0, 0)
    
    % Find the position of enzyme draws
    draw_pos = find( contains(model_ale.rxns, 'draw_'), 1, 'first');

    % Obtain flux fold changes and enzyme fold changes
    fluxfc = sol_ale.x(draw_pos:end-1) ./ sol_wt.x(draw_pos:end-1);
    exprfc = model_ale.enzFoldChange;

    % Remove NaN, Inf, zero fluxes enzymes.
    fluxfc(exprfc == 1) = [];
    exprfc(exprfc == 1) = [];
    exprfc(isnan(fluxfc) | (fluxfc == Inf)  | (fluxfc == 0)) = [];
    fluxfc(isnan(fluxfc) | (fluxfc == Inf)  | (fluxfc == 0)) = [];

    % Find the outliers
    if strcmpi(experiment,'nickel')
        perc.(experiment) = '100';
        fitted.(experiment) = fitlm(fluxfc, exprfc);
    else
    % Apply the generalized extreme Studentized deviate test for outliers
        outliers = find(isoutlier(fluxfc,'gesd'));
        perc.(experiment) = num2str(round(100 - length(outliers) / length(fluxfc) * 100, 2));
        fitted.(experiment) = fitlm(fluxfc, exprfc, 'Exclude', outliers);
    end
end 

%% ------------ FIGURES ------------
% First preview the figure, then set the paper size and color from the
% settings 'Print Preview'. Then call the print line.

%% Plot freq figure of the expressions
figure('Color',[1 1 1], 'WindowState', 'Maximized'); 
for exp = 1:length(experiments)
    experiment = experiments{exp} ;
    tt_exp = toptable.(experiment);
    subplot(2,4,exp); ggplot();
        histogram(tt_exp.logFC,'FaceAlpha', 0.8, 'FaceColor', ale_models.(experiment).color);
        n = length(tt_exp.logFC);
        ylabel(sprintf('Frequency (n=%d)',  n))
        xlabel('logFC') 
     title(ale_models.(experiment).legendName);
end
print('figures_pdf/expression_histograms_before','-dpdf', '-r600')
 
%% Plot freq figure of the expressions COMMONS
figure('Color',[1 1 1], 'WindowState', 'Maximized'); 
for exp = 1:length(experiments)
    experiment = experiments{exp} ;
    data = ale_models.(experiment).enzFoldChange;
    data(data == 1) = [];
    subplot(2,4,exp); ggplot();
    histogram(log2(data),'FaceAlpha', 0.8, 'FaceColor', ale_models.(experiment).color);
    n = ale_models.(experiment).enzFoldChangeN;
    title(ale_models.(experiment).legendName);
    ylabel(sprintf('Frequency (n=%d)',  n))
    xlabel('logFC')
end; clear n;
print('figures_pdf/expression_histograms_after','-dpdf', '-r600')

%%  Plot of the fitted data to see if the results are parallel
figure('Color',[1 1 1],'Position', [0 0 1280 1080]);
for exp = 1:length(experiments)    % Loop is blocked for the ethanol only
    subplot(2,4,exp); 
    
    experiment = experiments{exp} ;
    ggplot();
    fig = plot(fitted.(experiment)); 
    
    % First element is the data:
    set(fig(1), 'Color', [0.85,0.33,0.10], 'Marker', 'x', 'MarkerSize', 10);
    % Second element is the fitted line:
    set(fig(2), 'Color', 'black', 'LineWidth', 1, 'LineStyle', '--');
    % Third and forth elements are the confidence bounds:
    set(fig(3),'visible','off');
    set(fig(4),'visible','off');
    % Fifth is the x=y line
    line([0 2.5], [0 2.5], 'LineStyle', ':', 'Color', 'r') % 2.5 is selected for the ethanol, use max() for others if required
    
    % R text on the figure
    x = get(fig(2), 'XData'); 
    y = get(fig(2), 'YData'); 
    text(quantile(x, 0.9), quantile(y, 0.9), sprintf('R^{2} = %0.5f', fitted.(experiment).Rsquared.Ordinary), 'Interpreter', 'tex', 'FontSize', 14); 
    
    % Set title and labels
    title('')
    ylabel('Expression Fold Changes')
    xlabel('Enzyme-Draw Flux Fold Changes (v_{ale} / v_{wt})', 'Interpreter', 'tex');
    legend({sprintf('%s datapoints (%%%s)', ale_models.(experiment).legendName, perc.(experiment)), 'fitted line', 'x = y'}, 'Location','northoutside');
   
    hold off;
end 
print('figures_pdf/expression_fitted_lines','-dpdf', '-r600')
 