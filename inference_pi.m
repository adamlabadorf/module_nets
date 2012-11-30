function [pi_prim,pi_prim_samples] = inference_pi (options, assignment, program, data_binding)

% initialize pi_prim_params for sending to binding_likelihood
for ii = 1:options.num_modules
    pi_prim_params{ii} = rand(1,length(program{ii}.regulators));
end

max_ll = -10000;

for iter = 1:options.mh_samples
%for iter = 1:100

    for ii = 1:options.num_modules % run MH on each element of muprim sepeartely

        for jj = 1:length(program{ii}.regulators)

            % evaluate likelihood before proposal
            [orig_ll]= binding_likelihood(options,assignment,program,pi_prim_params,data_binding);

            % propose a new pi_prim for module ii and parent jj
            % pi_prim are probabilities, just draw from uniform (0,1)
            pi_prim_proposal = pi_prim_params;

            old_pi = pi_prim_params{ii}(jj);
            % I implemented the proposal pi as logistic_function(old_pi + N(0,1))
            % but it doesn't perform very well
            %new_pi = 1/(1+exp(-old_pi-randn(1)));

            % randomly +/-rand(0,1)*0.01 to old_pi
            new_pi = old_pi + (-1)^(binornd(1,0.5)+1)*rand(1)*0.01;
            new_pi = max(0,min(1,new_pi));

            pi_prim_proposal{ii}(jj) = new_pi;

            % evaluate likelihood after proposal
            [proposal_ll] = binding_likelihood(options,assignment,program,pi_prim_proposal,data_binding);

            lratio = exp(proposal_ll - orig_ll);
            accept_proposal = binornd(1,min(lratio,1));
            if accept_proposal == 1
                pi_prim_params{ii}(jj) = pi_prim_proposal{ii}(jj);
            end
        end
    end 

    pi_prims{iter}.pi_prim = pi_prim_params;

    [orig_ll]= binding_likelihood(options,assignment,program,pi_prim_params,data_binding);
    pi_prims{iter}.ll =  orig_ll;

    if orig_ll > max_ll
        max_ll = orig_ll;
        max_ll_iter = iter;
        pi_prim_max = pi_prim_params;
    end

end

pi_prim = pi_prim_max;
pi_prim_samples = pi_prims;
