function [options, modules, data_binding] = generate_data()

    options = create_options();
    modules = generate_modules(options);
    data_binding = generate_binding(options, modules);
