# This program run multiple WE simulation on multiple networks with different parameters

import os
import numpy as np

if __name__ == '__main__':

    # Netwrok parameters

    # N = [300,400,500,600,700,800,900,1000,1100,1200,1300,1400]
    N = 10000
    prog = 'gam'
    lam = 1.2
    # lam = 1+np.logspace(-2,0,9)
    # lam = np.array([1.5,1.6,1.7,1.8])
    measurements = 100000
    # eps_din = np.random.uniform(0.0, 3.0,measurements)
    # eps_din = [0.0, 0.05, 0.1, 0.15, 0.2]
    eps_din = 0.0
    eps_dout = eps_din
    # correlation = [-0.01,-0.03,-0.05,-0.08,-0.1,-0.12,-0.15,-0.18,-0.2,-0.25,-0.3]
    correlation = 0.0
    number_of_networks = 10
    # k = [50]

    k= 50

    # Catstrophe parameters parameters
    relaxation_time = 20
    sims = int(measurements/number_of_networks)
    tau = 100
    # tau = np.linspace(0.1,2.0,20)
    start = 50
    phi = 0.05
    # duartion = 1.0
    duartion = [0.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0]

    strength = 1.0-phi
    error_graphs = False

    # Parameters that don't change

    x = 0.2
    Alpha = 1.0

    # Paths needed to run the program
    dir_path = os.path.dirname(os.path.realpath(__file__))
    slurm_path = dir_path + '/slurm.serjob python3'
    program_path = dir_path + '/runwesim.py'
    loop_over = duartion
    run_mc_simulation = False

    for i in loop_over:
        error_graphs_flag = '--error_graphs' if error_graphs else ''
        run_mc_simulation_flag = '--run_mc_simulation' if run_mc_simulation else ''
        short_flag_flag = False
        command = (f'{slurm_path} {program_path} --N {N} --prog {prog} --lam {lam} --eps_din {eps_din} '
                   f'--eps_dout {eps_dout} --correlation {correlation} --number_of_networks {number_of_networks} '
                   f'--k {k} {error_graphs_flag} --sims {sims} --tau {tau} --start {start} --duartion {i} '
                   f'--strength {strength} --relaxation_time {relaxation_time} --x {x} '
                   f'--Alpha {Alpha} {run_mc_simulation_flag}')
        os.system(command)


        # for i in loop_over:
        #     error_graphs_flag = '--error_graphs' if error_graphs else ''
        #     run_mc_simulation_flag = '--run_mc_simulation' if run_mc_simulation else ''
        #     command = (
        #         f'{slurm_path} {program_path} --N {N} --prog {prog} --lam {lam} --eps_din {i} '
        #         f'--eps_dout {i} --correlation {correlation} --number_of_networks {number_of_networks} '
        #         f'--k {k} {error_graphs_flag} --sims {sims} --tau {tau} --start {start} --duartion {duartion} '
        #         f'--strength {strength} --relaxation_time {relaxation_time} --x {x} '
        #         f'--Alpha {Alpha} {run_mc_simulation_flag} {normalization_run_flag}'

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