function program = generate_program(options)

    regulators = options.regulators;

    for ii = 1:(options.num_modules)

        % each module has two regulators
        parents = regulators(ceil(rand(1,2)*length(regulators)));
        program{ii}.regulators = parents; 

    end

