function [modules] = infer_assignment_gibbs(options,modules,binding)

    model_ll = prob_bind(options,modules,binding);

    assignment = get_module_assignment(modules);

    for m = modules
        module_lls(m.id) = prob_bind_module(options,m,binding);
    end

    for gene_i = 1:options.num_genes
        orig_module = modules(assignment(gene_i));

        % get the module ll and compute its geometric mean
        orig_module_ll = module_lls(orig_module.id);
        mean_orig_ll = orig_module_ll/(length(orig_module.regulators)*length(orig_module.genes));

        % if gene is the only one in the module don't try to move it
        if length(orig_module.genes) == 1
            continue;
        end

        % find likelihood if this gene was removed from its module
        prop_module = orig_module;
        prop_module.genes(prop_module.genes==gene_i) = [];
        orig_gene_ll = prob_gene(gene_i,orig_module,binding);

        prop_ll = orig_module_ll - orig_gene_ll;
        mean_prop_ll = prop_ll/(length(prop_module.regulators)*length(prop_module.genes));

        % originating module has no change in probability
        ll_ratios(orig_module.id) = 0; 

        % for all other modules besides origin
        for m = setdiff(1:options.num_modules,orig_module.id)

            % get ll for this module
            other_ll = module_lls(other_module.id);

            % normalize the probability to the geometric mean
            mean_other_orig_ll = other_orig_ll/(length(modules(m).regulators)*length(modules(m).genes));

            % what would be the contribution of adding gene_i to this module?
            other_gene_ll = prob_gene(gene_i,other_module,binding);

            % concatenate the gene onto module
            other_module = modules(m);
            other_module.genes = [other_module.genes gene_i];
            other_modules(m) = other_module;

            other_prop_ll = other_ll + other_gene_ll;
            mean_other_prop_ll = other_orig_ll/(length(other_module.regulators)*length(other_module.genes));
            other_modules_lls(m_i) = other_prop_ll;

            ll_ratios(m_i) = (mean_prop_ll + mean_other_prop_ll) - (mean_orig_ll + mean_other_orig_ll);
            %ll_ratios(m) = other_gene_ll/length(other_module.regulators) - orig_gene_ll/length(orig_module.regulators);

        end

        % check case of a new module
        new_module = generate_module(options);
        %new_module.genes = [gene_i];

        % TODO - this needs to be adjusted for reversible jump or something
        %new_module_ll = prob_bind_module(new_module,binding);
        %ll_ratios(new_module.id) = (new_module_ll + ll_wo_gene) - (orig_module_ll + log(rand));

        % pick one of the modules in proportion to their probability ratios
        choice_proportions = exp(ll_ratios);
        choice = randsamp(1:length(ll_ratios),choice_proportions,1);
        %choice = [1:length(ll_ratios)](ll_ratios == max(ll_ratios))(1);
        if choice ~= orig_module.id % picked a different module
            %['moving ' num2str(gene_i) ' from ' num2str(orig_module.id) ' to ' num2str(choice)]

            % remove from old module
            modules(orig_module.id) = prop_module;
            module_lls(orig_module.id) = prop_module_ll;

            modules(choice) = other_modules(choice);
            module_lls(choice) = other_modules_lls(choice);

            % adjust model log likelihood
            model_ll += ll_ratios(choice);
        end
        assignment = get_module_assignment(modules);
        %sprintf('%d/%d\n',sum(assignment == real_assignment),length(assignment))
    end
end

function [gene_ll] = prob_gene(gene_i,module,binding)
        gene_binding = full(binding(module.regulators,gene_i));
        gene_ll = gene_binding.*log(module.pi_prim(module.regulators));
        gene_ll += (1-gene_binding).*log(1-module.pi_prim(module.regulators));
        gene_ll = sum(gene_ll);
end
