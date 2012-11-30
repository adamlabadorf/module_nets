function [binding,parameters]= generate_binding(options,assignment,program,pi_prim)

% INPUTS
%     assignment: assignment of genes into modules
%     program: regulatory program consisting of a number of regulators and a decision tree

if nargin < 4
    program = generate_program(options);
 
    if nargin < 3
        assignment = generate_assignment(options);
%         
%         if nargin < 2
%             parameters = options.parameters;
%         end
    end
end



for ii = 1:options.num_modules
    parents = program{ii}.regulators;
    
    for jj = 1:options.num_genes
        
        if assignment(jj) == ii % gene belongs to this module
            clear data;

            for kk = 1:length(parents)
                data(kk) = binornd(1,pi_prim{ii}(kk));
            end

            parameters{jj}.pi = pi_prim{ii};
            
            binding{jj}.data = data;
            binding{jj}.gene_id = jj;
            binding{jj}.module = ii;
            
            binding{jj}.num_modules = options.num_modules;
%             binding{jj}.assignment = assignment;
%             binding{jj}.program = program;
            
           % binding{jj}.parameters = parameters;
            
        end
    end
end
