function [modules] = infer_regulators_gibbs(options,modules,binding)

    model_ll = prob_bind(options,modules,binding);

    regulator_assignment = get_regulator_assignment(modules);

    for r = options.regulators
        orig_module = modules(regulator_assignment(r));
        reg_module_ll = prob_bind_module(options,orig_module,binding);

        % if removing regulator would destroy the module leave it
        if length(orig_module.regulators) == 1
            continue
        end

        prop_module = orig_module;
        prop_module.pi_prim(r) = 0;
        prop_module.regulators(prop_module.regulators == r) = [];
        %prop_module_ll = prob_bind_module(options,prop_module,binding);

        orig_reg_ll = prob_reg(r,orig_module,binding);

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

            %prop_module_ll = prob_bind_module(options,prop_module,binding);
            %other_ll = prob_bind_module(options,other_module,binding);

            other_reg_ll = prob_reg(r,other_module,binding);
            ll_ratios(m_i) = other_reg_ll/length(other_module.genes) - orig_reg_ll/length(orig_module.genes);

            %test_ratio = other_reg_ll - orig_reg_ll
            %real_ratio = (other_ll + prop_module_ll) - (other_orig_ll + reg_module_ll)

            %orig_ops = length(orig_module.regulators)*length(orig_module.genes);
            %prop_ops = length(prop_module.regulators)*length(prop_module.genes);
            %other_ops = length(modules(m_i).regulators)*length(modules(m_i).genes);
            %other_prop_ops = length(other_module.regulators)*length(other_module.genes);
            %diff_ops = (orig_ops-prop_ops)+(other_ops-other_prop_ops);

            % calculate log likelihood
            %prop_module_ll = prob_bind_module(options,prop_module,binding);
            %other_ll = prob_bind_module(options,other_module,binding);

            %orig_ratio = proposal_ratio(options,orig_module,modules(m_i),binding);
            %prop_ratio = proposal_ratio(options,prop_module,other_module,binding);
            %ll_ratios(m_i) = orig_ratio + prop_ratio;

            % log likelihood ratio over originating
            %[r orig_module.id m_i other_ll prop_module_ll other_orig_ll reg_module_ll];
            %ll_ratios(m_i) = (other_ll + prop_module_ll) - (other_orig_ll + reg_module_ll);

        end

        choice_proportions = exp(ll_ratios);
        choice = randsamp(1:length(ll_ratios),choice_proportions,1);
        %choice = [1:length(ll_ratios)](ll_ratios == max(ll_ratios))(1);
        if choice ~= orig_module.id % picked a different module
            %['moving ' num2str(gene_i) ' from ' num2str(gene_module.id) ' to ' num2str(choice)]

            % remove from old module
            modules(orig_module.id) = prop_module;

            if choice != orig_module.id % moved to existing
                modules(choice) = other_modules(choice);
            end

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
