
[options, program, real_assignment, real_pi_prim, data_binding] = generate_data();

% inference
[assignment] = inference_gene_assignment(options,assignment,program,data_binding);
