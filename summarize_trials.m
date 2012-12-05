#!/usr/bin/octave -qf
trial_fns = argv();

function [regs,genes] = compare_modules(mod1,mod2)
    for i = 1:length(mod1)
        for j = 1:length(mod2)
            regs(i,j) = length(mod1(i).regulators)-length(setdiff(mod1(i).regulators,mod2(j).regulators));
            genes(i,j) = length(mod1(i).genes)-length(setdiff(mod1(i).genes,mod2(j).genes));
        end
    end
end

for ii = 1:length(trial_fns)
    trial_fns{ii}
    load(trial_fns{ii});
    [regs, genes] = compare_modules(best_modules,real_modules)
end
