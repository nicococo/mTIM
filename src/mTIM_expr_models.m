function mTIM_expr_models(exp_name)
% Show the accuracy gain of a classifier during
% training (for paper: elegans_unfiltered)
% 

exp_name = 'elegans_unfiltered'
time_stamp = '130121_expr_models'

warning('This is a special function and not a general purpose tool.');

CFG = general_settings(exp_name,[],0,1);

sat = [1 2 3 5 10 20];
sat = [5 7 10 20];
sat = [15 20 5 3 2 1];
for i=1:length(sat),

    cls = sprintf('iter: %i',sat(i));
    CFG.PAR.num_levels = sat(i); 
    % Important!
    CFG.PAR.perm_feature_ranges = set_permitted_feature_ranges(CFG);
    
    main_run_training(CFG);
    CFG.eval_filename = sprintf('eval_expr_models%i',sat(i));
    main_run_prediction(CFG);
end

