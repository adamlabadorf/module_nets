function [modules] = infer_assignment_gibbs(options,modules,binding)

    model_ll = prob_bind(options,modules,binding);

    assignment = get_module_assignment(modules);

    for gene_i = 1:options.num_genes
        gene_module = modules(assignment(gene_i));

        % if gene is the only one in the module don't try to move it
        if length(gene_module.genes) == 1
            continue
        end

        % find likelihood if this gene was removed from its module
        module_wo_gene = gene_module;
        module_wo_gene.genes(module_wo_gene.genes==gene_i) = [];

        %ll_wo_gene = prob_bind_module(options,module_wo_gene,binding);
        %gene_ll = prob_gene(gene_i,gene_module,binding);

        % originating module has no change in probability
        ll_ratios(gene_module.id) = 0; 

        %gene_module_ll = prob_bind_module(options,gene_module,binding);
        orig_gene_ll = prob_gene(gene_i,gene_module,binding);
        %ll_wo_gene = gene_module_ll - orig_gene_ll;
        %ll_wo_gene_test = prob_bind_module(options,module_wo_gene,binding);
        % for all other modules besides origin
        for m = setdiff(1:options.num_modules,gene_module.id)
            %max_genes = max(length(gene_module.genes),length(other_module.genes));
            %gene_module_ll = prob_bind_module(options,gene_module,binding);
            %other_orig_ll = prob_bind_module(options,other_module,binding);

            % concatenate the gene onto module
            other_module = modules(m);
            other_module.genes = [other_module.genes gene_i];
            other_modules(m) = other_module;

            %other_ll = prob_bind_module(options,other_module,binding);
            other_gene_ll = prob_gene(gene_i,other_module,binding);
            %other_orig_ll = other_ll - other_gene_ll;
            %other_orig_ll_test = prob_bind_module(options,modules(m),binding);

            %orig_ratio = proposal_ratio(options,gene_module,modules(m),binding);
            %prop_ratio = proposal_ratio(options,module_wo_gene,other_module,binding);
            %ll_ratios(m) = orig_ratio + prop_ratio;
            ll_ratios(m) = other_gene_ll/length(other_module.regulators) - orig_gene_ll/length(gene_module.regulators);

            % calculate log likelihood
            %ll_wo_gene = prob_bind_module(options,module_wo_gene,binding);
            %other_ll = prob_bind_module(options,other_module,binding);

            % log likelihood ratio over originating
            %ll_ratios(m) = (other_ll + ll_wo_gene) - (other_orig_ll + gene_module_ll)
            %orig_ops = length(gene_module.regulators)*length(gene_module.genes);
            %prop_ops = length(module_wo_gene.regulators)*length(module_wo_gene.genes);
            %other_ops = length(modules(m).regulators)*length(modules(m).genes);
            %other_prop_ops = length(other_module.regulators)*length(other_module.genes);
            %diff_ops = (orig_ops-prop_ops)+(other_ops-other_prop_ops);

            %ll_ratios(m) = (other_gene_ll - orig_gene_ll) + diff_ops*log(0.5);

            % log likelhood ratio of just the gene reassignment contribution
            %ll_ratios(m) = (other_ll - other_orig_ll) - (gene_module_ll - ll_wo_gene);

            % log likelihood of the new model
            %ll_ratios(m) = (other_ll + ll_wo_gene);
        end

        % check case of a new module
        new_module = generate_module(options);
        %new_module.genes = [gene_i];

        % TODO - this needs to be adjusted for reversible jump or something
        %new_module_ll = prob_bind_module(new_module,binding);
        %ll_ratios(new_module.id) = (new_module_ll + ll_wo_gene) - (gene_module_ll + log(rand));

        % pick one of the modules in proportion to their probability ratios
        choice_proportions = exp(ll_ratios);
        choice = randsamp(1:length(ll_ratios),choice_proportions,1);
        %choice = [1:length(ll_ratios)](ll_ratios == max(ll_ratios))(1);
        if choice ~= gene_module.id % picked a different module
            %['moving ' num2str(gene_i) ' from ' num2str(gene_module.id) ' to ' num2str(choice)]

            % remove from old module
            modules(gene_module.id) = module_wo_gene;

            if choice == new_module.id % new module
                modules(new_module.id) = new_module;
                options.num_modules += 1;
            elseif choice != gene_module.id % moved to existing
                modules(choice) = other_modules(choice);
            end

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
