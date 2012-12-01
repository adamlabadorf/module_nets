function [options, program, assignment, pi_prim, data_binding] = generate_data()

    options = create_options();
    program = generate_program(options);
    assignment = generate_assignment(options);
    pi_prim = generate_pi_prim(options,program);
    data_binding = generate_binding(options, assignment, program, pi_prim);
