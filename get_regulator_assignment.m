function [assignment] = get_regulator_assignment(modules)
    for m = modules
        assignment(m.regulators) = m.id;
    end
end
