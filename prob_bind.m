function [ll] = prob_bind(options,modules,binding)
    for module = modules
        ll_m(module.id) = prob_bind_module(options,module,binding);
    end
    ll = sum(ll_m);
end


