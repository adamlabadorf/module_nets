function [ options ] = create_options()

options.num_conditions = 10;
options.num_genes = 100;
options.num_regulators = 10;
options.num_modules = 5;
options.regulators = 1:10;

options.expression_model = 'normal';
parameters.mu_prim = 0;
parameters.sigma_prim = 1;
parameters.pi_prim = 0.1;
 
options.parameters = parameters;
options.sigma_coeff = 1;

options.mh_samples = 1000;
end

