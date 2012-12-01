function [assignment] = inference_gene_assignment (options,assignment,program,pi_prim_params,data_binding)

[orig_ll]= binding_likelihood(options,assignment,program,pi_prim_params,data_binding);

    % Initial State
    for iter = 1:options.mh_samples
        for ii = 1:options.num_genes

            old_num_modules = options_old.num_modules;

            %%%%% Metropolis-Hastings %%%%%

            % picke the module to move
            old_module = assignment_old(ii);

            %%% Proposal:

            % randomly assign the gene to another module or to its own new
            % cluster

            proposal_module = ceil(rand(1,1)*old_num_modules);
            proposal_num_modules = old_num_modules;

            % % %     % if the random is to its own module, use it as case of new cluster
            % % %     %%% have to change this, should use reverse jump mcmc for this
            % % %     if proposal_module == old_module
            % % %         proposal_module = old_num_modules + 1;
            % % %         proposal_num_modules = proposal_module; %(last index=number)
            % % %     end
            % % %     
            assignment_proposal = assignment_old;
            assignment_proposal(ii) = proposal_module;

            %%% Likelihood ratio  
            % Should change this to only calculate the likelihood change for the
            % target module, not to recalculate the overall likelihood

            exp_param_old = exp_param;
            exp_param_old.assignment = assignment_old;
            likelihood_old = (model_likelihood(options_old,exp_param,bind_param,data_expression,data_binding));
            exp_param_proposal = exp_param;
            exp_param_proposal.assignment = assignment_proposal;
            options_proposal = options_old;
            options_proposal.num_modules = proposal_num_modules;
            likelihood_proposal = (model_likelihood(options_proposal,exp_param_proposal,bind_param,data_expression,data_binding));

            lratio = exp(likelihood_proposal - likelihood_old);

            % We assume the proposal is symmetric so Q(muprim_old|muprim_prim) = Q(muprim_prim|muprim_old)
            % is this TRUE in this case?

            accept_proposal = binornd(1,min(lratio,1));
            if accept_proposal == 1
                assignment_new = assignment_proposal;
            end

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

        %%%%% check back this function, might have issues with proper updating of
        %%%%% the values and retrning them back to the param structs
