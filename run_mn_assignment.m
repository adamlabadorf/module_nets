
[options, real_modules, binding] = generate_data();

% reassign regulators for inference
trial_modules = real_modules;
for mm = 1:options.num_modules
    trial_modules(mm).regulators = options.regulators(ceil(rand(1,2)*length(options.regulators)));
end

% inference
[modules] = inference_regulator_assignment(options, trial_modules, binding);
