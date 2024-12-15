# This program run multiple WE simulation on multiple networks with different parameters

import os
import numpy as np

if __name__ == '__main__':

    # Netwrok parameters

    # N = [300,400,500,600,700,800,900,1000,1100,1200,1300,1400]
    N = 1200
    prog = 'bd'
    lam = 1.18
    # lam = 1+np.logspace(-2,0,9)
    # lam = np.array([1.5,1.6,1.7,1.8])
    measurements = 100000
    # eps_din = np.random.uniform(0.0, 3.0,measurements)
    eps_din = [0.0, 0.05, 0.1, 0.15, 0.2]
    eps_dout = eps_din
    # correlation = [-0.01,-0.03,-0.05,-0.08,-0.1,-0.12,-0.15,-0.18,-0.2,-0.25,-0.3]
    correlation = 0.0
    number_of_networks = 100
    # k = [50]

    k= 80

    # Catstrophe parameters parameters
    relaxation_time = 20
    sims = int(measurements/number_of_networks)
    tau = 150
    # tau = np.linspace(0.1,2.0,20)
    start = 50
    phi = 1.0
    duartion = 1.0
    strength = 1.0-phi
    error_graphs = False

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
        short_flag_flag = False
        command = (f'{slurm_path} {program_path} --N {N} --prog {prog} --lam {lam} --eps_din {i} '
                   f'--eps_dout {i} --correlation {correlation} --number_of_networks {number_of_networks} '
                   f'--k {k} {error_graphs_flag} --sims {sims} --tau {tau} --start {start} --duartion {duartion} '
                   f'--strength {strength} --relaxation_time {relaxation_time} --x {x} '
                   f'--Alpha {Alpha} {run_mc_simulation_flag}')
        os.system(command)