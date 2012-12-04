function [modules] = infer_regulators_gibbs(options,modules,binding)

    model_ll = prob_bind(options,modules,binding);

    regulator_assignment = get_regulator_assignment(modules);

    for m = modules
        module_lls(m.id) = prob_bind_module(options,m,binding);
    end

    for r = options.regulators
        orig_module = modules(regulator_assignment(r));
        %orig_ll = prob_bind_module(options,orig_module,binding);
        orig_ll = module_lls(orig_module.id);
        mean_orig_ll = orig_ll/(length(orig_module.regulators)*length(orig_module.genes));

        % if removing regulator would destroy the module leave it
        if length(orig_module.regulators) == 1
            continue
        end

        prop_module = orig_module;
        prop_module.pi_prim(r) = 0;
        prop_module.regulators(prop_module.regulators == r) = [];
        %prop_ll = prob_bind_module(options,prop_module,binding);

        orig_reg_ll = prob_reg(r,orig_module,binding);
        prop_ll = orig_ll - orig_reg_ll;
        mean_prop_ll = prop_ll/(length(prop_module.regulators)*length(prop_module.genes));

        % no change to likelhood if we stay the same
        %ll_ratios(orig_module.id) = 0;
        for m_i = setdiff(1:options.num_modules,orig_module.id)

            %gene_module_ll = prob_bind_module(options,orig_module,binding);
            %other_orig_ll = prob_bind_module(options,modules(m_i),binding);

            % concatenate the gene onto module
            other_module = modules(m_i);
            other_module.pi_prim(r) = orig_module.pi_prim(r);
            other_module.regulators = [other_module.regulators r];
            other_modules(m_i) = other_module;

            %prop_ll = prob_bind_module(options,prop_module,binding);
            %other_ll = prob_bind_module(options,other_module,binding);

            other_ll = module_lls(other_module.id);
            mean_other_ll = other_ll/(length(modules(m_i).regulators)*length(modules(m_i).genes));

            other_reg_ll = prob_reg(r,other_module,binding);
            other_prop_ll = other_ll + other_reg_ll;
            mean_other_prop_ll = other_prop_ll/(length(other_module.regulators)*length(other_module.genes));
            other_modules_lls(m_i) = other_prop_ll;

            %ll_ratios(m_i) = other_reg_ll/length(other_module.genes) - orig_reg_ll/length(orig_module.genes);
            ll_ratios(m_i) = (mean_prop_ll + mean_other_prop_ll) - (mean_orig_ll + mean_other_ll);

        end

        choice_proportions = exp(ll_ratios);
        choice = randsamp(1:length(ll_ratios),choice_proportions,1);
        %choice = [1:length(ll_ratios)](ll_ratios == max(ll_ratios))(1);
        if choice ~= orig_module.id % picked a different module
            %['moving ' num2str(gene_i) ' from ' num2str(gene_module.id) ' to ' num2str(choice)]

            % remove from old module
            modules(orig_module.id) = prop_module;
            module_lls(orig_module.id) = prop_ll;

            % moved to existing
            modules(choice) = other_modules(choice);
            module_lls(choice) = other_modules_lls(choice);

            % adjust model log likelihood
            model_ll += ll_ratios(choice);
        end
    end
end

function [reg_ll] = prob_reg(reg,module,binding)
        reg_binding = full(binding(reg,module.genes));
        reg_ll = reg_binding.*log(module.pi_prim(reg));
        reg_ll += (1-reg_binding).*log(1-module.pi_prim(reg));
        reg_ll = sum(reg_ll);
end
