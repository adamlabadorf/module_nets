function [ratio] = proposal_ratio(options,orig,prop,binding)

    ll_orig = prob_bind_module(options,orig,binding);
    ll_prop = prob_bind_module(options,prop,binding);

    orig_ops = length(orig.regulators)*length(orig.genes);
    prop_ops = length(prop.regulators)*length(prop.genes);

    ll_orig /= orig_ops;
    ll_prop /= prop_ops;

    %diff_ops = orig_ops - prop_ops;
    %if sign(diff_ops) == -1 % fewer opts in orig, adjust
    %    ll_orig += abs(diff_ops)*log(0.5);
    %elseif sign(diff_ops) == 1
    %    ll_prop += diff_ops*log(0.5);
    %end

    ratio = ll_prop - ll_orig;

end
