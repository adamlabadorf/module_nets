#!/usr/bin/octave -qf
out_fn = argv(){1};
kind = argv(){2};
% assignment
% 
% Definitions:
%     N - # genes
%     M - # modules
%     R - # regulators
% 
% assign genes randomly to modules
% assign regulators randomly to modules
% 
% P(B|S,A) = \prod_{m \in M} \prod_{r \in m} \prod_{g \in G_m} (# bound)/(total)
% 
% while not converged:
% 
%     1. foreach g of N genes:
%            calculate probability of model as is
%            calculate probability of model without g
%            foreach m of M modules:
%                calculate probability of model when reassigning gene g to module m
%            make a module with p random regulators
%            calculate probability of assigning g to random module as only gene
%            randomly select one of the M+1 modules in proportion to the probability ratios
%            if there is only one gene left in a module, leave it as is
%              - later might want to actually reassign it with reversible jump
%     2. foreach r of R regulators:
%            calculate probability of model as is
%            foreach m of M modules :
%                if r \not \in m :
%                    calculate probability of model when r is assigned to module m
%                else :
%                    calculate probability of model when r is removed from m
%            randomly select one of the M modules in proportion to the probability ratios
%            add/remove regulator from module
%            if all regulators moved out of module:
%                redistribute genes randomly to other modules

function [modules] = randomize_assignment(options,in_modules)
    rand_modules = generate_modules(options);
    for ii = 1:length(in_modules)
        modules(ii) = in_modules(ii);
        modules(ii).genes = rand_modules(ii).genes;
    end
end

function [modules] = randomize_pi(options,modules)
    for m_i = 1:length(modules)
        modules(m_i).pi_prim = zeros(length(modules(m_i).regulators),1);
        modules(m_i).pi_prim(modules(m_i).regulators) = 0.5+rand(length(modules(m_i).regulators),1)*0.5;
    end
end

function [modules] = randomize_regulators(options,modules)
    regulators = options.regulators;
    % assign regulators by sampling w/o replacement
    for mm = 1:length(modules)
        new_module = modules(mm);
        new_module.regulators = [];
        new_module.pi_prim = zeros(length(options.regulators),1); % make all pi_prims 0
        for ii = 1:options.init_modules
            if length(regulators) != 0
                new_reg = randsamp(regulators,ones(length(regulators),1),1);
                new_module.regulators(ii) = new_reg;
                new_module.pi_prim(new_reg) = modules(mm).pi_prim(modules(mm).regulators(ii));
                regulators = setdiff(regulators,new_module.regulators);
            end
        end
        modules(mm) = new_module;
    end
end 

function [regs,genes] = compare_modules(mod1,mod2)
    for i = 1:length(mod1)
        for j = i:length(mod2)
            regs(i,j) = length(setdiff(mod1(i).regulators,mod2(j).regulators));
            genes(i,j) = length(setdiff(mod1(i).genes,mod2(j).genes));
        end
    end
end

options = create_options(kind);

[binding,real_modules] = generate_binding(options);
real_ll = prob_bind(options,real_modules,binding)
real_assignment = get_module_assignment(real_modules);

% randomize assignment
modules = real_modules;
modules = randomize_assignment(options,modules);

% randomize pi
%modules = randomize_pi(options,modules);

% randomize regulators
modules = randomize_regulators(options,modules);

assignment = get_module_assignment(modules);

before_module = modules;
before_ll = prob_bind(options,modules,binding)

[regs,genes] = compare_modules(before_module,real_modules)

for iter = 1:options.mh_samples

    % infer assignment
    modules = infer_assignment_gibbs(options,modules,binding);
    num_genes = arrayfun(@(x) length(x.genes),modules);

    % infer pi
    %modules = infer_pi_mh(options,modules,binding);
    %modules = infer_pi_direct(options,modules,binding);

    % infer regulators
    modules = infer_regulators_gibbs(options,modules,binding);

    ll = prob_bind(options,modules,binding);

    lls(iter) = ll;
    module_samples{iter} = modules;

end


after_ll = prob_bind(options,modules,binding)

best_iteration = [1:length(lls)](lls == max(lls));
best_modules = module_samples{lls == max(lls)};
[regs,genes] = compare_modules(best_modules,real_modules)
'real regulators'
real_modules.regulators
'best regulators'
best_modules.regulators
'real genes'
real_modules.genes
'best genes'
best_modules.genes


save(argv(){1})
