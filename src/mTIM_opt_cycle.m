function mTIM_opt_cycle(exp_name,dir_name)
% Show the accuracy gain of a classifier during
% training (for paper: elegans_unfiltered)
% 

warning('This is a special function and not a general purpose tool.');


sat = [2 5 10 15 20 30 45 60 80 100 200 500];
%sat = [25 35 50 70];
%sat = [80 85];
for i=1:length(sat),

    cls = sprintf('iter: %i',sat(i));
    CFG = general_settings(exp_name,dir_name,0,0);

    % set the classifier iter 
    CFG.iteration = sat(i);
    CFG.eval_filename = sprintf('eval_iter%i',sat(i));
 
    main_run_prediction(CFG);
end

