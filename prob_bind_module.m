function [ll] = prob_bind_module(options,module,binding)
    for r = module.regulators
        % the genes that are bound by r
        binding_events = full(binding(r,module.genes));

        % p(module|binding) = p(binding|module)p(module)/p(binding)
        p = module.pi_prim(r);
        num_bound = sum(binding_events);
        num_unbound = sum(1-binding_events);
        ll(r) = num_bound*log(p)+num_unbound*log(1-p);

        % calculate the integral for p(binding)
        f = @(x) x^num_bound * (1-x)^num_unbound;
        ll(r) -= log(quad(f,0,1));
    end

    ll = sum(ll);

end
