function [ options ] = create_options(kind)

if nargin == 0
    kind = 'easy'
end

if strcmp(kind,'easy')
    options.num_genes = 100;
    options.num_modules = 5;
    options.init_modules = 2;

    options.regulator_names = 1:10;
    options.num_regulators = length(options.regulator_names);
    options.regulators = 1:options.num_regulators;
    options.bind_prob = 0.99;
    options.mh_samples = 1000;

elseif strcmp(kind,'medium')
    options.num_genes = 300;
    options.num_modules = 10;
    options.init_modules = 3;

    options.regulator_names = 1:30;
    options.num_regulators = length(options.regulator_names);
    options.regulators = 1:options.num_regulators;
    options.bind_prob = 0.85;
    options.mh_samples = 5000;

elseif strcmp(kind,'hard')
    options.num_genes = 1000;
    options.num_modules = 25;
    options.init_modules = 4;

    options.regulator_names = 1:100;
    options.num_regulators = length(options.regulator_names);
    options.regulators = 1:options.num_regulators;
    options.bind_prob = 0.8;
    options.mh_samples = 3;
end

