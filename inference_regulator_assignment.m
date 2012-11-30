function [program] = inference_regulator_assignment(options, assignment, program, pi_prim, data_binding)

% check out structure learning in bayes nets for this

all_regulators = options.regulators;
options_old = options;


% Real parent assignment
program_real = exp_param.program;

% Initialize randomly
program_old = generate_program(options);


for iter = 1:options.mh_samples
    for mm = 1:options.num_modules
        
        % real paren assignment for that module
        module_parents_real = program{mm}.regulators;

        %%% Proposal:
        % randomly add or remove a regulator to/from the parents of the module
        % if chosen to add a new regulaotr, we should check the condition for
        % acyclity

        add_regulator = binornd(1,0.5);

        if add_regulator == 1
            not_in_parents = setdiff(all_regulators,module_parents_real);
            to_add = randsample(not_in_parents,1);
            module_parents_proposal = [module_parents_old, to_add];
        else
            remove_ind = ceil(rand(1,1)*length(module_parents_real);
            module_parents_proposal = [module_parents_old(1:remove_ind-1), module_parents_old(remove_ind+1:end)];
        end

        program_proposal = program_old;
        assignment_proposal(ii) = proposal_module;

        %%% Likelihood ratio  
        % Should change this to only check the likelihood of the specific
        % module, no tto re-calculate the whole

        lratio = exp(likelihood_proposal - likelihood_old);


        % We assume the proposal is symmetric so Q(muprim_old|muprim_prim) = Q(muprim_prim|muprim_old)
        % is this TRUE in this case?
        
        if lratio>=1
            assignment_new = assignment_proposal;

        else
            accept_proposal = binornd(1,lratio);
            if accept_proposal == 1
                assignment_new = assignment_proposal;
            else
                assignment_new = assignment_old;
            end
        end


    %     assignment_samples(:,ii + (iter-1)*length(muprim_old)) = assignment_new;
    %     ratio(ii,iter) = lratio;

        assignment_old = assignment_new;

        samples(:,ii + (iter-1)*options.num_genes) = assignment_new;
        
    end
            
    exp_param.assignment = assignment_old;
    data_expression.parameters = exp_param;

    options = options_old;

end

figure,
for ii = 1:5
    subplot(5,1,ii), hist(samples(ii,:),50);
    hold on, plot(assignment_real(ii),0,'*r','MarkerSize',10), 
end

       
      
data_expression.parameters = exp_param;
data_binding.parameters = bind_param;


        
        
        
        
        
        
        
        

        
