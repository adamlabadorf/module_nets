function [binding,modules] = generate_binding(options,modules)

% binding is a PxG binary matrix
% - P is # of all regulators
% - G is # of genes

if nargin == 1
    modules = generate_modules(options);
end

binding = sparse(options.num_regulators,options.num_genes);

for module = modules
    parents = module.regulators;
    %pi_prim = module.pi_prim;
    for gene = module.genes
        for r = options.regulators
            if ismember(r,parents) == 1
                binding(r,gene) = binornd(1,options.bind_prob);
            else
                binding(r,gene) = binornd(1,1-options.bind_prob);
            end
        end
    end
end
