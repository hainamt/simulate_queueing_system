% 5 simulations
sim_indices = [1 2 3 4 5];
total_customer = [17 12 11 19 6]';
K = 15;
entering_sims = sim_indices(total_customer(sim_indices) < K);
server_states = [0 0; % sim 1 free sv1 and sv2
                 1 1;
                 0 1; % sim 3 free sv1
                 1 0; % sim 4 free sv2
                 2 1];

is_available = server_states == 0;

