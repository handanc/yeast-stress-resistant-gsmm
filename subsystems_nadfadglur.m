clear;clc;
load ale_models
experiments = { 'ethanolb2','ethanolb8', 'caffeine', 'coniferylaldehyde', 'iron', 'nickel', 'phenylethanol', 'silver', 'wildtype'};
figlegend = {'B2 (ethanol)', 'B8 (ethanol)', 'CAF905-2 (caffeine)', 'BH-13 (coniferyl aldehyde)', 'M8FE (iron)', 'M9 (nickel)', 'C9 (phenylethanol)', '2E (silver)', 'Reference'}

load fba_results.mat
load fva_results.mat


%% Define Subsystems
namelist = { 'sce00010', % Glycolysis
            'sce00020', % Citrate cycle (TCA cycle)
            'sce00190', % Oxidative phosphorylation
            'sce00030', % Pentose phosphate pathway
           % 'sce01200', % Carbon metabolism
            'sce00620', % Pyruvate metabolism
            'sce01210', % 2-Oxocarboxylic acid metabolism
            'sce00600', % Sphingolipid metabolism
           % 'sce00630', % Glyoxylate and dicarboxylate metabolism
            'sce01212', % Fatty acid metabolism
            %'sce00230', % Purine metabolism
            %'sce00240', % Pyrimidine metabolism
            'sce01230', % Biosynthesis of amino acids
            'sce01110', % Biosynthesis of secondary metabolites
            %'sce01130', % Biosynthesis of antibiotics
        }
    
namelistSubs = {'Glycolysis', 
            'Citrate cycle (TCA cycle)',
            'Oxidative phosphorylation',
            'Pentose phosphate pathway',
            %'Carbon metabolism',
            'Pyruvate metabolism',
            '2-Oxocarboxylic acid metabolism',
            'Sphingolipid metabolism',
           % 'Glyoxylate and dicarboxylate metabolism',
            'Fatty acid metabolism',
            %'Purine metabolism',
            %'Pyrimidine metabolism',
            'Biosynthesis of amino acids',
            'Biosynthesis of secondary metabolites',
            %'Biosynthesis of antibiotics'
        }

% Find rxns for specific subsystems
[subSystems, subsys] = listSubsystems(ale_models.wildtype, namelist);

FVAsub = [];
for i=1:length(namelist)
    sub = namelist{i};
    FVAsub.span.(sub).all = [];
    for exp = 1:length(experiments)
            experiment = experiments{exp} ;
            FVAsub.span.(sub).(experiment) = FVA.(experiment)(subsys.(sub),3) - FVA.(experiment)(subsys.(sub),2);
            FVAsub.span.(sub).all = [FVAsub.span.(sub).all;  FVAsub.span.(sub).(experiment)];
            FVAsub.avg.(sub)(1,exp) = mean(FVAsub.span.(sub).(experiment)(FVAsub.span.(sub).(experiment) < 1000));
            FVAsub.avg.all(i,exp) = mean(FVAsub.span.(sub).(experiment)(FVAsub.span.(sub).(experiment) < 1000)); 
    end
end

%% Plot: Vertical Stacked bar chart

sums = sum(FVAsub.avg.all);
percentages = (FVAsub.avg.all ./ sums) *100;

figure('Color',[1 1 1]); ggplot();
colors = jet(length(namelistSubs));
%colors = { '#4e46f6','#1a8dec', '#37ca4e', '#fcba03', '#fc4948', '#fe7e0e', '#f429c1', '#8240e3'};
%colors = hex2rgb(colors);
percentages = percentages';
H = bar(percentages,'stacked'); 
ylabel('Subsystem Flux Span Percentage')
xticks(1:length(experiments)+1)
xtickangle(45)
xticklabels(experiments)
for i = 1:length(colors)
   H(i).FaceColor = 'flat';
    H(i).CData = colors(i,:);
end
 
% A loop that does num2str conversion only if value is > 1%
for i=1:size(percentages,1);
    for j=1:size(percentages,2);
        if percentages(i,j)>2;
        labels_stacked=num2str(percentages(i,j),'%.1f %%');
        hText = text(i, sum(percentages(i,1:j), 2), labels_stacked);
        set(hText, 'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize', 12, 'Color','k');
        end
    end
end
ylim([0 100])
a = legend(namelistSubs,'Location', 'eastoutside', 'box', 'off')

%% Plot: Horizontal Stacked bar chart

colors = { '#885953','#e15a58', '#ee9a52', '#ffd740', '#97cc8f', '#e8e8e8', '#3f9ed8', '#b2d3ee', '#ccaed7', '#9f68b3'};
colors = hex2rgb(colors);

percentages = percentages'; %

% whole section below must be run at once
figure('Color',[1 1 1], 'WindowState', 'maximized'); 
subplot(3,1,1:2); ggplot();
    grid off;
    H = barh(percentages, 'stacked'); 
    yticks(1:length(experiments)+1)
    yticklabels(figlegend)
    xlim([0 100])
    xticks([0:5:100])
    legend(namelistSubs,'Location', 'southoutside', 'box', 'off', 'NumColumns', 5)

    for i = 1:length(colors)
        H(i).FaceColor = 'flat';    % sometime throws error, rerun fixes, idk why.
        H(i).CData = colors(i,:);
    end
    set(gca, 'XAxisLocation','top')
    a = get(gca,'XTickLabel');
    set(gca,'fontsize',15)

% couldn't find a way to automatically add numbers onto plot so writing
% numbers onto bars from below table is done on the adobe acrobat.
tableplot = subplot(3,1,3);
tableplot.XColor = 'white';
tableplot.YColor = 'white';
uitable('Data', percentages, 'Position', [500 120 800 200]);

% Save the figure
print('figures_pdf/f8_a_subsytems','-dpdf', '-r600')


%% Print normalized numbers as heatmap
normalized = (FVAsub.avg.all ./ FVAsub.avg.all(:,9)); 
h = heatmap(normalized(:,1:8))
h.XDisplayLabels = figlegend(1:8);
h.YDisplayLabels = namelistSubs; 
print('figures_pdf/f8_b_subsystem_heatmap','-dpdf', '-r600')


%% Collect fluxes in a single vector
model = ale_models.wildtype;
rxnFluxes = [];
for exp = 1:length(experiments)
    experiment = experiments{exp} ; 
    rxnFluxes(:,exp) = fba_solvectors.gur10.(experiment);
end

%% ATP Production and Consumption reactions
metATP = {'s_0434'};
metATP = {'s_0434', 's_0435', 's_0437', 's_0438', 's_0439'};

atpProRxns = findRxnsFromMets(model, metATP, false, 'producersOnly', true);
atpConRxns = findRxnsFromMets(model, metATP, false, 'consumersOnly', true);

atpProRxns = findRxnIDs(model, atpProRxns);
atpConRxns = findRxnIDs(model, atpConRxns);
 

% Producing
rxnFluxes_Pro = rxnFluxes(atpProRxns,:);
rxnStds_Pro = std(rxnFluxes_Pro, [], 2);
rxnList_Pro = model.rxnNames(atpProRxns);
rxnTable_Pro = table(atpProRxns, rxnList_Pro, rxnFluxes_Pro, rxnStds_Pro);
rxnTable_Pro(all(rxnTable_Pro.rxnFluxes_Pro == 0, 2),:) = [];
rxnTable_Pro = sortrows(rxnTable_Pro, 'rxnStds_Pro', 'descend');

% Consuming
rxnFluxes_Con = rxnFluxes(atpConRxns,:);
rxnStds_Con = std(rxnFluxes_Con, [], 2);
rxnList_Con = model.rxnNames(atpConRxns);
rxnTable_Con = table(atpConRxns, rxnList_Con, rxnFluxes_Con, rxnStds_Con);
rxnTable_Con(all(rxnTable_Con.rxnFluxes_Con == 0, 2),:) = [];
rxnTable_Con = sortrows(rxnTable_Con, 'rxnStds_Con', 'descend');


sums.gur10.atp_pro = sum(rxnTable_Pro.rxnFluxes_Pro);
sums.gur10.atp_con = sum(rxnTable_Con.rxnFluxes_Con); 

%% NAD reactions

metsNAD = { 's_1198', 's_1199', 's_1200', 's_1201', 's_1202', ...
            's_1203', 's_1204', 's_1205', 's_1206', 's_1207', ...
            's_1208', 's_1210', 's_1211', 's_1212', 's_1213', ...
            's_1214', 's_1215'};
nadProRxns = findRxnsFromMets(model, metsNAD, false, 'producersOnly', true);
nadConRxns = findRxnsFromMets(model, metsNAD, false, 'consumersOnly', true);

nadProRxns = findRxnIDs(model, nadProRxns);
nadConRxns = findRxnIDs(model, nadConRxns);

% Producing
rxnFluxes_Pro = rxnFluxes(nadProRxns,:);
rxnStds_Pro = std(rxnFluxes_Pro, [], 2);
rxnList_Pro = model.rxnNames(nadProRxns);
rxnTable_Pro = table(nadProRxns, rxnFluxes_Pro, rxnFluxes_Pro, rxnStds_Pro);
rxnTable_Pro(all(rxnTable_Pro.rxnFluxes_Pro == 0, 2),:) = [];
rxnTable_Pro = sortrows(rxnTable_Pro, 'rxnStds_Pro', 'descend');


% Consuming
rxnFluxes_Con = rxnFluxes(nadConRxns,:);
rxnStds_Con = std(rxnFluxes_Con, [], 2);
rxnList_Con = model.rxnNames(nadConRxns);
rxnTable_Con = table(nadConRxns, rxnFluxes_Con, rxnFluxes_Con, rxnStds_Con);
rxnTable_Con(all(rxnTable_Con.rxnFluxes_Con == 0, 2),:) = [];
rxnTable_Con = sortrows(rxnTable_Con, 'rxnStds_Con', 'descend');

sums.gur10.nad_pro = sum(rxnTable_Pro.rxnFluxes_Pro);
sums.gur10.nad_con = sum(rxnTable_Con.rxnFluxes_Con);

%% FAD
metsFAD = { 's_0687', 's_0688', 's_0689', 's_0690' };

fadProRxns = findRxnsFromMets(model, metsFAD, false, 'producersOnly', true);
fadConRxns = findRxnsFromMets(model, metsFAD, false, 'consumersOnly', true);

fadProRxns = findRxnIDs(model, fadProRxns);
fadConRxns = findRxnIDs(model, fadConRxns);


rxnFluxes_Pro = rxnFluxes(fadProRxns,:);
rxnStds_Pro = std(rxnFluxes_Pro, [], 2);
rxnList_Pro = model.rxnNames(fadProRxns);
rxnTable_Pro = table(fadProRxns, rxnList_Pro, rxnFluxes_Pro, rxnStds_Pro);
rxnTable_Pro(all(rxnTable_Pro.rxnFluxes_Pro == 0, 2),:) = [];
rxnTable_Pro = sortrows(rxnTable_Pro, 'rxnStds_Pro', 'descend');


rxnFluxes_Con = rxnFluxes(fadConRxns,:);
rxnStds_Con = std(rxnFluxes_Con, [], 2);
rxnList_Con = model.rxnNames(fadConRxns);
rxnTable_Con = table(fadConRxns, rxnList_Con, rxnFluxes_Con, rxnStds_Con);
rxnTable_Con(all(rxnTable_Con.rxnFluxes_Con == 0, 2),:) = [];
rxnTable_Con = sortrows(rxnTable_Con, 'rxnStds_Con', 'descend');

sums.gur10.fad_pro = rxnTable_Pro.rxnFluxes_Pro;
sums.gur10.fad_con = rxnTable_Con.rxnFluxes_Con;


%% Glutathione metabolism

[gsa, gsaPos]  = findRxnsFromSubSystem(model,{'sce00480  Glutathione metabolism'});
 
rxnFluxes_Glu = rxnFluxes(gsaPos,:);
rxnStds_Glu = std(rxnFluxes_Glu, [], 2);
rxnList_Glu = model.rxnNames(gsaPos);
rxnTable_Glu = table(gsaPos, rxnList_Glu, rxnFluxes_Glu, rxnStds_Glu);
rxnTable_Glu(all(rxnTable_Glu.rxnFluxes_Glu == 0, 2),:) = [];
rxnTable_Glu = sortrows(rxnTable_Glu, 'rxnStds_Glu', 'descend');

sums.gur10.glut = sum(rxnTable_Glu.rxnFluxes_Glu);

%% Save Table
save('data/fba_sumtables.mat', 'sums')

%% NAD FAD GLUT TABLE

freeplot = [sums.gur10.nad_con; sums.gur10.fad_con; sums.gur10.glut]';

figure('Color',[1 1 1], 'WindowState','maximized');ggplot()
bar(freeplot, 'BarWidth', 1, 'EdgeColor', 'none')
xticks(1:9)
xticklabels(multiLineLabels())
ylabel('Total Flux (mmol/gDWh)')
legend('NAD, NADH, NADP, NADPH incl. reactions', 'FAD, FADH2 incl. reactions', 'Glutathione Metabolism reactions')

print('figures_pdf/f7_nad_fad_glut','-dpdf', '-r600')

%% Amino acid

[~, aminoacidsPos]  = findRxnsFromSubSystem(model,{'sce01230  Biosynthesis of amino acids'});
 
rxnFluxes_AA = rxnFluxes(aminoacidsPos,:);
rxnStds_AA = std(rxnFluxes_AA, [], 2);
rxnList_AA = model.rxnNames(aminoacidsPos);
rxnTable_AA = table(aminoacidsPos, rxnList_AA, rxnFluxes_AA, rxnStds_AA);
rxnTable_AA(all(rxnTable_AA.rxnFluxes_AA == 0, 2),:) = [];
rxnTable_AA(contains(rxnTable_AA.rxnList_AA, '(No'),:) = [];
rxnTable_AA(~contains(rxnTable_AA.rxnList_AA, 'syn'),:) = [];
rxnTable_AA = sortrows(rxnTable_AA, 'rxnStds_AA', 'descend');

%% 12 precursor metabolites

metspre = { 's_0568', ... %d-glucose-6-phosphate (G6P), 
            's_0557', ... %d-fructose-6-phosphate (F6P), 
            's_0577', ... %d-ribulose-5-phosphate (R5P), 
            's_0551', ... %d-erythrose-4-phosphate (E4P), 
            's_0764', ... %glyceraldehyde-3-phosphate (GAP),
            's_0260', ... %3-phosphonato-D-glycerate (3PG), 
            's_1360', ... %phosphoenolpyruvate (PEP),
            's_1399', ... %pyruvate (PYR), 
            's_0373', ... %acetyl-CoA (ACA), 
            's_0180', ... %2-oxoglutarate (2KG)
            's_1464', ... %succinyl-CoA (SCA)
            's_1271', ... %oxaloacetate (OXA)
            };

%% LOOPS REQUIRED BELOW    <<<<<<<<<<<<<<   UNFINISHED 

metspre = {'s_1271'};

prePro = findRxnsFromMets(model, metspre, false, 'producersOnly', true);
preCon = findRxnsFromMets(model, metspre, false, 'consumersOnly', true);

prePro = findRxnIDs(model, prePro);
preCon = findRxnIDs(model, preCon);

prepre = preCon;

rxnList = model.rxnNames(prepre);
rxnEnz = findGenesFromRxns(model, model.rxns(prepre));
rxnFluxes = vectors.gur10(prepre,2:end);
rxnStds = std(rxnFluxes, [], 2);
rxnTable_Pro = table(prepre, rxnList,rxnEnz, rxnFluxes, rxnStds);
rxnTable_Pro(all(rxnTable_Pro.rxnFluxes == 0, 2),:) = [];
rxnTable_Pro = sortrows(rxnTable_Pro, 'rxnStds', 'descend');


% FBA for biomass precursors
surfNet(model, model.rxns(754))
% 1229, 'r_4041', biomass pseudoreaction; 754 r_2111 growth
metabolist = { 's_1096',   %lipid
                's_3717',   %protein
                's_3718',   %carbohydrate
                's_3719',   %RNA
                's_3720'}   %DNA

metabolistID = findMetIDs(model, metabolist)
model.metNames(metabolistID)


precursors = findRxnsFromMets(model, model.mets(metabolistID), false, 'producersOnly', true);
%    {'r_2108'}
%     {'r_4047'}
%     {'r_4048'}
%     {'r_4049'}
%     {'r_4050'}
    
%Protein pseudoreaction
aminoacids = findMetsFromRxns(model, 'r_4047')
aminoacidsPro = findRxnsFromMets(model, aminoacids, false, 'producersOnly', true);
aminoacidsProID = findRxnIDs(model,aminoacidsPro);



for exp = 1:length(experiments) 
    experiment = experiments{exp} ;
    for i=1:length(aminoacidsProID)
        tmp = ale_models.(experiment);
        tmp = changeObjective(tmp, tmp.rxns(aminoacidsProID(i)));
        sol = optimizeCbModel(tmp, 'max'); 
        vecf.(experiment)(i,1) = sol.f;
        vecf.(experiment)(i,2) = sol.x(754);
    end
end

for exp = 1:length(experiments) 
    experiment = experiments{exp} ;
    for i=1:length(aminoacidsProID)
        tmp = ale_models.(experiment);
        sol = optimizeCbModel(tmp, 'max');
        veckof.(experiment)(i,1) = sol.f;
        tmp = changeRxnBounds(tmp, tmp.rxns(aminoacidsProID(i)), 0 , 'b');
        sol = optimizeCbModel(tmp, 'max'); 
        veckof.(experiment)(i,2) = sol.f;
    end
end



%Lipid pseudoreaction
lipids_backbone = findMetsFromRxns(model, 'r_4063');
lipids_chain = findMetsFromRxns(model, 'r_4065');
lipids_backbonePro = findRxnsFromMets(model, lipids_backbone, false, 'producersOnly', true);
lipids_backboneProID = findRxnIDs(model,lipids_backbonePro);
lipids_chainPro = findRxnsFromMets(model, lipids_chain, false, 'producersOnly', true);
lipids_chainProID = findRxnIDs(model,lipids_chainPro);

% lipid chain
for exp = 1:length(experiments) 
    experiment = experiments{exp} ;
    for i=1:length(lipids_backboneProID)
        tmp = ale_models.(experiment);
        tmp = changeObjective(tmp, tmp.rxns(lipids_backboneProID(i)));
        sol = optimizeCbModel(tmp, 'max'); 
        vecf.(experiment)(i,1) = sol.f;
        if any(sol.x)
        vecf.(experiment)(i,2) = sol.x(754);
        end
    end
end

for exp = 1:length(experiments) 
    experiment = experiments{exp} ;
    for i=1:length(lipids_backboneProID)
        tmp = ale_models.(experiment);
        sol = optimizeCbModel(tmp, 'max');
        veckof.(experiment)(i,1) = sol.f;
        tmp = changeRxnBounds(tmp, tmp.rxns(lipids_backboneProID(i)), 0 , 'b');
        sol = optimizeCbModel(tmp, 'max'); 
        veckof.(experiment)(i,2) = sol.f;
    end
end

