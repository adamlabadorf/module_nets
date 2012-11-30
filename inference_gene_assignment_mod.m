function [ll]= binding_likelihood(options,assignment,program,pi_prim,B,module)

% parameters = B.parameteres;
% pi_prim = parameters.pi_prim;

pi_module = pi_prim{module};

for jj = 1:options.num_genes
    
    if assignment(jj) == ii % gene belongs to this module
        
       gene_data = B{jj}.data;
       ll_gene(jj) = sum(gene_data .* log(pi_module) + (1-gene_data).*log(1-pi_module));
       
    end
end

ll = sum(ll_gene);
