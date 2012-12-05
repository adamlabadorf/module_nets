function [modules] = infer_regulators_gibbs(options,modules,binding)

    model_ll = prob_bind(options,modules,binding);

    regulator_assignment = get_regulator_assignment(modules);

    for m = modules
        module_lls(m.id) = prob_bind_module(options,m,binding);
    end

    regs_to_update = 1:options.num_regulators;
    %regs_to_update = randsamp(1:options.num_regulators,ones(options.num_regulators),floor(options.num_regulators*0.1));
    for r =regs_to_update 
        orig_module = modules(regulator_assignment(r));
        %orig_ll = prob_bind_module(options,orig_module,binding);
        %orig_ll = module_lls(orig_module.id);
        orig_ll = prob_bind_module(options,orig_module,binding);

        % if removing regulator would destroy the module leave it
        if length(orig_module.regulators) == 1
            continue
        end

        prop_module = orig_module;
        prop_module.pi_prim(r) = 0;
        prop_module.regulators(prop_module.regulators == r) = [];
        %prop_ll = prob_bind_module(options,prop_module,binding);

        %orig_reg_ll = prob_reg(r,orig_module,binding);
        %prop_ll = orig_ll - orig_reg_ll;
        prop_ll = prob_bind_module(options,prop_module,binding);

        % no change to likelhood if we stay the same
        ll_ratios(orig_module.id) = 0;
        ll(orig_module.id) = orig_ll;
        for m_i = setdiff(1:options.num_modules,orig_module.id)

            %other_orig_ll = module_lls(other_module.id);
            other_orig_ll = prob_bind_module(options,modules(m_i),binding);

            % concatenate the gene onto module
            other_module = modules(m_i);
            other_module.pi_prim(r) = orig_module.pi_prim(r);
            other_module.regulators = [other_module.regulators r];
            other_modules(m_i) = other_module;

            %prop_ll = prob_bind_module(options,prop_module,binding);
            other_prop_ll = prob_bind_module(options,other_module,binding);

            %other_reg_ll = prob_reg(r,other_module,binding);
            %other_prop_ll = other_orig_ll + other_reg_ll;
            other_prop_ll = prob_bind_module(options,other_module,binding);
            other_modules_lls(m_i) = other_prop_ll;

            ll_ratios(m_i) = (prop_ll + other_prop_ll) - (orig_ll + other_orig_ll);
        end

        choice_proportions = exp(ll_ratios);
        choice = randsamp(1:length(choice_proportions),choice_proportions,1);
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
