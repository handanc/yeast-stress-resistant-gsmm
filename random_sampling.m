%% Load Model Files
clear;clc;
load ale_models
experiments = { 'ethanolb2','ethanolb8', 'caffeine', 'coniferylaldehyde', 'iron', 'nickel', 'phenylethanol', 'silver', 'wildtype'};

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


%% Sampling
for exp = 1:length(experiments)
    
    experiment = experiments{exp};
    tmp = ale_models.(experiment);
    %tmp.ub(tmp.ub == Inf) = 1000;
    tmp = changeRxnBounds(tmp, exchange_rxns_uptakes, 0, 'b');
    tmp = changeRxnBounds(tmp, exchange_rxns_uptakes_allowed, Inf, 'u'); 
    
    sol = optimizeCbModel(tmp, 'max', 0, 0);
    growth = sol.f;

    tmp = changeRxnBounds(tmp, tmp.rxns(753), sol.f, 'u'); 
    tmp = changeRxnBounds(tmp, tmp.rxns(753), sol.f*0.9, 'l'); 
    
    [solutions90.(experiment) goodRxns.(experiment)] = randomSampling(tmp, 5000, true, false, true );
end
save('data/randomsampling_90th_5k', 'solutions90', 'goodRxns')

points = solutions90;