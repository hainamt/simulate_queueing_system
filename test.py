import numpy as np

# Small example demonstrating handle_a2_events vectorization
def handle_a2_example():
    print("=== Handle A2 Events Vectorization Example ===\n")

    # Setup: 5 simulations, 2 servers, capacity K=4
    num_sims = 5
    C = 2  # number of servers
    K = 4  # capacity

    # Initialize system state
    total_customers = np.array([1, 3, 4, 2, 0])  # Current customers in each sim
    server_states = np.array([
        [1, 2],  # sim 0: server 0 busy (S1), server 1 busy (S2)
        [0, 1],  # sim 1: server 0 idle, server 1 busy (S1)
        [2, 2],  # sim 2: both servers busy (S2)
        [0, 0],  # sim 3: both servers idle
        [0, 0]   # sim 4: both servers idle
    ])

    # Simulations having A2 events (indices 0, 2, 4)
    sim_indices = np.array([0, 2, 4])

    print("Initial State:")
    print(f"sim_indices with A2 events: {sim_indices}")
    print(f"total_customers: {total_customers}")
    print(f"server_states:\n{server_states}")
    print(f"K (capacity): {K}")
    print()

    # Step 1: Check which customers can enter
    print("Step 1: Check which can enter the system")
    can_enter = total_customers[sim_indices] < K
    print(f"total_customers[sim_indices]: {total_customers[sim_indices]}")
    print(f"can_enter mask: {can_enter}")

    entering_sims = sim_indices[can_enter]
    print(f"entering_sims: {entering_sims}")
    print()

    # Step 2: Update customer count for entering customers
    print("Step 2: Update customer count")
    total_customers[entering_sims] += 1
    print(f"total_customers after increment: {total_customers}")
    print()

    # Step 3: Vectorized server assignment
    print("Step 3: Vectorized server assignment")

    # Extract server states for entering simulations
    server_states_entering = server_states[entering_sims]
    print(f"server_states_entering:\n{server_states_entering}")

    # Find idle servers (state == 0)
    idle_mask = server_states_entering == 0
    print(f"idle_mask (True = idle server):\n{idle_mask}")

    # Check which simulations have at least one idle server
    has_idle_server = np.any(idle_mask, axis=1)
    print(f"has_idle_server: {has_idle_server}")

    # Find the first idle server for each simulation
    first_idle_indices = np.argmax(idle_mask, axis=1)
    print(f"first_idle_indices (argmax of idle_mask): {first_idle_indices}")

    # Assign servers only where available
    server_assignments = np.full(len(entering_sims), -1, dtype=int)
    server_assignments[has_idle_server] = first_idle_indices[has_idle_server]
    print(f"server_assignments (-1 = no idle server): {server_assignments}")
    print()

    # Step 4: Update server states for valid assignments
    print("Step 4: Update server states")
    valid_mask = server_assignments >= 0
    valid_entering_sims = entering_sims[valid_mask]
    valid_server_assignments = server_assignments[valid_mask]

    print(f"valid_mask: {valid_mask}")
    print(f"valid_entering_sims: {valid_entering_sims}")
    print(f"valid_server_assignments: {valid_server_assignments}")

    if len(valid_entering_sims) > 0:
        print(f"\nBefore server state update:")
        print(f"server_states:\n{server_states}")

        # Update server states vectorially
        server_states[valid_entering_sims, valid_server_assignments] = 1

        print(f"After server state update (assigned servers set to 1):")
        print(f"server_states:\n{server_states}")

        # Calculate S1 event indices for scheduling
        s1_event_indices = 2 + 2 * valid_server_assignments
        print(f"s1_event_indices for scheduling: {s1_event_indices}")

    print()

    # Step 5: Handle losses (simulations that couldn't enter)
    print("Step 5: Handle losses")
    blocked = total_customers[sim_indices] >= K
    blocked_sims = sim_indices[blocked]
    print(f"total_customers[sim_indices] after processing: {total_customers[sim_indices]}")
    print(f"blocked mask (>= K): {blocked}")
    print(f"blocked_sims: {blocked_sims}")

    print()
    print("=== Final State ===")
    print(f"total_customers: {total_customers}")
    print(f"server_states:\n{server_states}")
    print(f"Simulations that got service: {valid_entering_sims}")
    print(f"Simulations that were blocked: {blocked_sims}")


handle_a2_example()