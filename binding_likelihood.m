function [ll2]= binding_likelihood(options,assignment,program,pi_prim,B)

% parameters = B.parameteres;
% pi_prim = parameters.pi_prim;

for ii = 1:options.num_modules
    
    
    pi_module = pi_prim{ii};
    
    for jj = 1:options.num_genes
        
        if assignment(jj) == ii % gene belongs to this module
            
           gene_data = B{jj}.data;
           ll_gene(jj) = sum(gene_data .* log(pi_module) + (1-gene_data).*log(1-pi_module));
           
        end
    end
end

ll2 = sum(ll_gene);
