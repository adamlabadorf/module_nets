function [module] = generate_module(options)
    module.id = options.num_modules+1;
    module.regulators = options.regulators(ceil(rand(options.init_modules,1)*length(options.regulators)));
    module.pi_prim = [];
end
