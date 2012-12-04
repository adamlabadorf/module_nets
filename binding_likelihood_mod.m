function [ll2]= binding_likelihood_mod(options,binding,module)

pi_module = module.pi_prim;

for gene_id = module.assignment
    
   gene_data = binding(module.regulators,gene_id)';
   ll_gene(gene_id) = sum(gene_data .* log(pi_module) + (1-gene_data).*log(1-pi_module));
   
end

[ll2] = sum(ll_gene);
