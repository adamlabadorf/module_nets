function pi_prim = generate_pi_prim(options,program);

for ii = 1:options.num_modules
    num_parents = length(program{ii}.regulators);
    pi_prim{ii} = rand(1,num_parents);
end
