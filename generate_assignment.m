function [assignment] = generate_assignment(options)

assignment = ceil(rand(options.num_genes,1)*options.num_modules);