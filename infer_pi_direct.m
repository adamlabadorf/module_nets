function [modules] = infer_pi_direct(options,modules,binding)
    for m_i = 1:length(modules)
        m = modules(m_i);
        for r = m.regulators
            % we put the min(0.99,x) in there because probability==1 causes numerical problems
            modules(m_i).pi_prim(r) = min(0.99,sum(binding(r,m.genes))/length(m.genes));
        end
    end
end
