function [modules] = infer_assignment_gibbs(options,modules,binding)

    model_ll = prob_bind(options,modules,binding);

    assignment = get_module_assignment(modules);

    for m = modules
        module_lls(m.id) = prob_bind_module(options,m,binding);
    end

    % how about just optimizing a random subset?
    genes_to_update = 1:options.num_genes;
    %genes_to_update = randsamp(1:options.num_genes,ones(options.num_genes),floor(options.num_genes*0.1));
    %genes_to_update = [1];
    for gene_i = genes_to_update
        orig_module = modules(assignment(gene_i));

        % get the module ll and compute its geometric mean
        %orig_ll = module_lls(orig_module.id);
        orig_ll = prob_bind_module(options,orig_module,binding);

        % if gene is the only one in the module don't try to move it
        if length(orig_module.genes) == 1
            continue;
        end

        % find likelihood if this gene was removed from its module
        prop_module = orig_module;
        prop_module.genes(prop_module.genes==gene_i) = [];
        orig_gene_ll = prob_gene(gene_i,orig_module,binding);
        prop_ll = orig_ll - orig_gene_ll;
        prop_ll = prob_bind_module(options,prop_module,binding);

        % originating module has no change in probability
        ll_ratios(orig_module.id) = 0; 

        % for all other modules besides origin
        for m_i = setdiff(1:options.num_modules,orig_module.id)

            other_module = modules(m_i);

            % get ll for this module
            %other_orig_ll = module_lls(other_module.id);
            other_orig_ll = prob_bind_module(options,other_module,binding);

            % what would be the contribution of adding gene_i to this module?
            %other_gene_ll = prob_gene(gene_i,other_module,binding);

            % concatenate the gene onto module
            other_module.genes = [other_module.genes gene_i];
            other_modules(m_i) = other_module;

            %other_prop_ll = other_orig_ll + other_gene_ll;
            other_prop_ll = prob_bind_module(options,other_module,binding);
            other_modules_lls(m_i) = other_prop_ll;

            ll_ratios(m_i) = (prop_ll + other_prop_ll) - (orig_ll + other_orig_ll);

        end

        % pick one of the modules in proportion to their probability ratios
        choice_proportions = exp(ll_ratios);
        choice = randsamp(1:length(choice_proportions),choice_proportions,1);
        if choice ~= orig_module.id % picked a different module
            %['moving ' num2str(gene_i) ' from ' num2str(orig_module.id) ' to ' num2str(choice)]

            % remove from old module
            modules(orig_module.id) = prop_module;
            module_lls(orig_module.id) = prop_ll;

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

function [frac] = frac_bind(gene_i,module,binding)
        gene_binding = full(binding(module.regulators,gene_i));
end
