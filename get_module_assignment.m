function [assignment] = get_module_assignment(modules)
    for module = modules
        assignment(module.genes) = module.id;
    end
end

