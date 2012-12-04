function [modules] = inference_regulator_assignment(options, modules, data_binding)

% check out structure learning in bayes nets for this

all_regulators = options.regulators;
options_old = options;


% we will use ll to adjust the overall likelihood by adding/subtracting each
% module likelihood
[ll] = binding_likelihood(options, data_binding, modules);

ll

for iter = 1:options.mh_samples
    for module = modules 
        
        proposal_module = module;

        %%% Proposal:
        % randomly add or remove a regulator to/from the parents of the module
        % if chosen to add a new regulaotr, we should check the condition for
        % acyclity

        add_regulator = binornd(1,0.5);

        if add_regulator == 1
            not_in_parents = setdiff(all_regulators,module.regulators);
            to_add = not_in_parents(ceil(rand*length(not_in_parents)));
            proposal_module.regulators = [module.regulators, to_add];
            proposal_module.pi_prim = [module.pi_prim rand]; % just add a random pi?
            % in this case, should we infer a good value of the new pi before calculating
            % the overall likelihood?
        else
            r = rand(1,1)
            remove_ind = ceil(r*length(module.regulators));
            proposal_module.regulators = module.regulators([1:remove_ind-1, remove_ind+1:length(module.regulators)]);
            proposal_module.pi_prim = module.pi_prim([1:remove_ind-1, remove_ind+1:length(module.regulators)]);
        end

        %%% Likelihood ratio
        ll_old = binding_likelihood_mod(options, data_binding, module);
        ll_proposal = binding_likelihood_mod(options, data_binding, proposal_module);

        % adjust for reversible jump
        if add_regulator == 1 % real is on bottom
            u = randn;
            p_u = normpdf(u);
            dg = @(x) exp(-x)./(1+exp(-x)).^2;
            lratio = ll_proposal + log(dg(u)) - ll_old - log(p_u);
        else
            remove_ind
            module.pi_prim(remove_ind)
            u = logit(module.pi_prim(remove_ind));
            p_u = normpdf(u);
            dg = @(x) -1/(x^2-x)
            lratio = ll_proposal + log(dg(u)) - ll_old + log(p_u);
        end

        accept = binornd(1,min(1,exp(lratio)));
        if accept == 1
            modules(module.id) = proposal_module;
        end

    end % end module loop

end % end MH loop
