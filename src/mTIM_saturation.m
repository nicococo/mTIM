function mTIM_saturation(exp_name)

warning('This tool is currently under development (do not use it yet)');


mydate = datestr(now,'yymmdd')
sat = [80 120 250 500 750 1000 1500]

for i=1:length(sat),

    dir_name = sprintf('%s_%i',mydate,sat(i));
    CFG = general_settings(exp_name,dir_name,0,1);

    % set the training example size 
    for j=1:size(CFG.train_params,1),
        CFG.train_params{j,4} = sat(i);
    end;

    if (i~=1),
        main_run_training(CFG);
    end
    CFG = general_settings(exp_name,dir_name,0,0);
    main_run_prediction(CFG);
end

