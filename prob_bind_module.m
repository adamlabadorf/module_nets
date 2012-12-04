function [ll] = prob_bind_module(options,module,binding)
    ll = 0;
    for r = module.regulators
        % the genes that are bound by r
        binding_events = full(binding(r,:))(module.genes);

        % ratio of bound/total
        %if sum(binding_events) > 0
        %    ll += log(sum(binding_events)) - log(length(module.genes));
        %end

        % genes that bind
        ll += sum(binding_events) * log(module.pi_prim(r));
        % genes that don't bind
        ll += sum(1-binding_events) * log(1-module.pi_prim(r));

    end

end
