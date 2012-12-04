function [ll2]= binding_likelihood(options,binding,modules)

for module = modules

    ll_module(module.id) = binding_likelihood_mod(options,binding,module);

end

ll2 = sum(ll_module);
