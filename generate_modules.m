function [modules] = generate_modules(options)

    assignment = generate_assignment(options);
    regulators = options.regulators;

    for mm = 1:options.num_modules

        modules(mm).id = mm;

        % generate program
        % each module starts with exactly two random parents
        %modules(mm).regulators = options.regulators(ceil(rand(options.init_modules,1)*length(options.regulators)));
        % assign regulators by sampling w/o replacement
        for ii = 1:options.init_modules
            if length(regulators) != 0
                modules(mm).regulators(ii) = randsamp(regulators,ones(length(regulators),1),1);
                regulators = setdiff(regulators,modules(mm).regulators);
            end
        end

        % generate pi_prim
        % uniform random of length of module's regulators
        modules(mm).pi_prim = zeros(length(options.regulators),1);
        %modules(mm).pi_prim(modules(mm).regulators) = options.bind_prob;
        modules(mm).pi_prim(modules(mm).regulators) = options.bind_prob + rand(length(modules(mm).regulators),1)*(1-options.bind_prob);

        modules(mm).genes = [1:options.num_genes](assignment == mm);

    end
end
