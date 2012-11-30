function program = generate_program(options)
    regulators = options.regulators;

for ii = 1:(options.num_modules)
    
    %parents = find_parents(module);
    parents = regulators(ceil(rand(1,2)*length(regulators)));
  
   % num_leaves = length(parents)+1;
    num_conditions = options.num_conditions;
    
%     leaf_conditions = partition_conditions(num_leaves,num_conditions);
%     
%     for part = 1:length(num_leaves)
%         
%     end

    %program{ii}.condition_state = (sign(rand(length(parents),num_conditions)-0.5)+1)/2;
    %program{ii}.condition_state = (sign(rand(length(parents),num_conditions)-0.5));
    program{ii}.regulators = parents; 
end
  
