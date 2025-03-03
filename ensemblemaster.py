# This program run multiple WE simulation on multiple networks with different parameters

import os
import numpy as np
import time

if __name__ == '__main__':

    # Netwrok parameters

    # N = [300,400,500,600,700,800,900,1000,1100,1200,1300,1400]
    N = 1000
    prog = 'gam'
    lam = 1.24
    # lam = 1+np.logspace(-2,0,9)
    # lam = np.array([1.5,1.6,1.7,1.8])
    # eps_din = np.random.uniform(0.0, 3.0,measurements)
    # eps_din = [0.0, 0.05, 0.1, 0.15, 0.2]
    eps_din = np.linspace(0.01,1.0,20)
    eps_dout = eps_din
    # measurements = np.where(eps_din < 0.2, 1000000, 10000)
    measurements = 1000000
    # correlation = [-0.01,-0.03,-0.05,-0.08,-0.1,-0.12,-0.15,-0.18,-0.2,-0.25,-0.3]
    correlation = 0.0
    number_of_networks = 50
    # k = [50]
    k= 50

    # Catstrophe parameters parameters
    relaxation_time = 20
    # sims = (measurements/number_of_networks).astype(int)
    sims = int(measurements/number_of_networks)
    tau = 100
    # tau = np.linspace(0.1,2.0,20)
    start = 50
    # phi = np.linspace(0.01,1.0,2)
    phi = 1.0
    duartion = 2.0
    # duartion = [0.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0]
    # duartion = np.linspace(0.0,5.0,20)

    strength = 1.0-phi
    # strength = np.ones(len(phi)) - phi
    error_graphs = False

    # Set this flag as needed:
    normalization_run = True  # Set to True if you want normalization to run

    # Create the flag string based on the value of normalization_run
    normalization_run_flag = '--normalization_run' if normalization_run else ''

    # Parameters that don't change

    x = 0.2
    Alpha = 1.0

    # Paths needed to run the program
    dir_path = os.path.dirname(os.path.realpath(__file__))
    slurm_path = dir_path + '/slurm.serjob python3'
    program_path = dir_path + '/runwesim.py'
    loop_over = eps_din
    run_mc_simulation = False

    for i in loop_over:
        error_graphs_flag = '--error_graphs' if error_graphs else ''
        run_mc_simulation_flag = '--run_mc_simulation' if run_mc_simulation else ''
        command = (
            f'{slurm_path} {program_path} --N {N} --prog {prog} --lam {lam} --eps_din {i} '
            f'--eps_dout {i} --correlation {correlation} --number_of_networks {number_of_networks} '
            f'--k {k} {error_graphs_flag} --sims {sims} --tau {tau} --start {start} --duartion {duartion} '
            f'--strength {strength} --relaxation_time {relaxation_time} --x {x} '
            f'--Alpha {Alpha} {run_mc_simulation_flag} {normalization_run_flag}'
        )
        # Initialize variables before first attempt
        result = 1  # Initialize to failure state
        retry_count = 0
        backoff_time = 2  # Start with a short delay

        while result != 0:
            if retry_count > 0:
                print(f"Retry #{retry_count} for eps_din={i}, waiting {backoff_time} seconds...")
                time.sleep(backoff_time)
                # Exponential backoff - increase wait time with each failure
                backoff_time = min(backoff_time * 2, 30)  # Cap at 30 seconds

            result = os.system(command)
            retry_count += 1

            if retry_count > 20:  # Prevent infinite loops
                print(f"Failed to submit job for eps_din={i} after 20 attempts. Skipping.")
                break

    # for d in duartion:
    #     for i, j in zip(loop_over, sims):
    #         error_graphs_flag = '--error_graphs' if error_graphs else ''
    #         run_mc_simulation_flag = '--run_mc_simulation' if run_mc_simulation else ''
    #         command = (
    #             f'{slurm_path} {program_path} --N {N} --prog {prog} --lam {lam} --eps_din {i} '
    #             f'--eps_dout {i} --correlation {correlation} --number_of_networks {number_of_networks} '
    #             f'--k {k} {error_graphs_flag} --sims {j} --tau {tau} --start {start} --duartion {d} '
    #             f'--strength {strength} --relaxation_time {relaxation_time} --x {x} '
    #             f'--Alpha {Alpha} {run_mc_simulation_flag} {normalization_run_flag}'
    #         )
    #         os.system(command)