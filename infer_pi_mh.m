function [modules] = infer_pi_mh(options,modules,binding)
    for iter = 1:options.mh_samples
        sprintf('infer_pi_mh: %d',iter);
        model_ll = prob_bind(modules,binding)
        for m_i = 1:length(modules)
            m = modules(m_i);
            for r = m.regulators

                orig_ll = prob_bind_module(m,binding);

                proposal_m = m;

                % randomly +/-rand(0,1)*0.01 to old_pi
                new_pi = m.pi_prim(r) + (-1)^(binornd(1,0.5)+1)*rand(1)*0.01;
                new_pi = max(0,min(1,new_pi));
                proposal_m.pi_prim(r) =  new_pi;

                proposal_ll = prob_bind_module(proposal_m,binding);

                lratio = exp(proposal_ll - orig_ll);
                accept_proposal = binornd(1,min(lratio,1));
                if accept_proposal == 1
                    modules(m_i) = proposal_m;
                end
            end
        end
    end
end
