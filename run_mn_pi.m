
[options, program, assignment, real_pi_prim, data_binding] = generate_data();

[likelihood] = binding_likelihood(options,assignment,program,real_pi_prim,data_binding)

% inference
[pi_prim_hat,pi_prim_samples] = inference_pi(options,assignment,program,data_binding);

pi_prim_hat

range = 1:length(pi_prim_samples);
for iter = range
    residuals(iter) = 0;
    for ii = 1:options.num_modules
        residuals(iter) += sum((real_pi_prim{ii}-pi_prim_samples{iter}.pi_prim{ii}).^2);
    end
    ll(iter) = pi_prim_samples{iter}.ll;
end

figure,
[AX,H1,H2] = plotyy(range,residuals,range,ll);
set(get(AX(1),'Ylabel'),'String','sum(residuals)^2');
set(get(AX(2),'Ylabel'),'String','binding Log-likelihood');
title('pi prime');
