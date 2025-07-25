# This program run multiple WE simulation on multiple networks with different parameters

import os
import numpy as np
import time
import rand_networks
import runwesim
import networkx as nx
import pickle
from scipy.stats import skew
import argparse
import time
from scipy.sparse.linalg import eigsh


if __name__ == '__main__':

    # Netwrok parameters

    # N = [300,400,500,600,700,800,900,1000,1100,1200,1300,1400]
    N = 1000
    prog = 'bd'
    lam = 1.6
    # lam = 1+np.logspace(-2,0,9)
    # lam = np.array([1.5,1.6,1.7,1.8])
    # eps_din = np.random.uniform(0.0, 3.0,measurements)
    # eps_din = [0.0, 0.05, 0.1, 0.15, 0.2]
    # eps_din = np.linspace(0.01, 1.0, 5)
    eps_din = 0.5
    eps_dout = eps_din
    # measurements = 1000000
    # correlation = [-0.01,-0.03,-0.05,-0.08,-0.1,-0.12,-0.15,-0.18,-0.2,-0.25,-0.3]
    # correlation = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
    # correlation = 0.0
    # correlation = np.linspace(-0.6, 0.6, 2)
    correlation = np.linspace(0.01, 0.5, 10)
    number_of_networks = 100
    # k = [50]
    k= 50

    # Catstrophe parameters parameters
    relaxation_time = 20
    # sims = (measurements/number_of_networks).astype(int)
    # sims = int(measurements/number_of_networks)
    tau = 100
    # tau = np.linspace(0.1,2.0,20)
    start = 50
    # phi = np.linspace(0.01,1.0,2)
    phi = 1.0
    duartion = 2.0
    # duartion = [0.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0]
    # duartion = np.linspace(0.01,10.0,20)

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
    run_mc_simulation = False
    heatmap = False
    correlation_heat_map = False

    # if correlation_heat_map:
    #     G = rand_networks.configuration_model_undirected_graph_mulit_type(k, eps_din, N, prog,0.0)

    def submit_job(N, prog, lam, eps_din, eps_dout, correlation, number_of_networks, k,
                   error_graphs, sims, tau, start, duartion, strength, relaxation_time, x,
                   Alpha, run_mc_simulation, normalization_run_flag,
                   slurm_path='',
                   program_path='',
                   runheatcorrelation=True,
                   infile='GNull.pickle'):
        error_graphs_flag = '--error_graphs' if error_graphs else ''
        run_mc_simulation_flag = '--run_mc_simulation' if run_mc_simulation else ''
        attempts = 20

        k = int(k)
        tau = int(tau)
        start = int(start)
        sims = int(sims)

        error_graphs_flag = '--error_graphs' if error_graphs else ''
        run_mc_simulation_flag = '--run_mc_simulation' if run_mc_simulation else ''
        runheat_flag = '--runheatcorrelation' if runheatcorrelation else ''

        command = (
            f"{slurm_path} {program_path}"
            f" --N {N} --prog {prog} --lam {lam}"
            f" --eps_din {eps_din} --eps_dout {eps_dout}"
            f" --correlation {correlation}"
            f" --number_of_networks {number_of_networks} --k {k}"
            f" {error_graphs_flag}"
            f" --sims {sims} --tau {tau} --start {start}"
            f" --duartion {duartion} --strength {strength}"
            f" --relaxation_time {relaxation_time}"
            f" --x {x} --Alpha {Alpha}"
            f" {run_mc_simulation_flag}"
            f" {normalization_run_flag}"
            f" {runheat_flag}"
            f" --graphname {infile}"
        )

        result = 1
        retry_count = 0
        backoff_time = 0.1

        while result != 0:
            if retry_count > 0:
                print(f"Retry #{retry_count} for eps_din={eps_din}, waiting {backoff_time} seconds...")
                time.sleep(backoff_time)
                backoff_time = min(backoff_time * 2, 0.5)

            result = os.system(command)
            retry_count += 1

            if retry_count > attempts:
                print(f"Failed to submit job for eps_din={eps_din} after {attempts} attempts. Skipping.")
                break

    def submit_correlation_heatmap(G,correlation_heat_map,parameters):

        def submit_with_retries(slurm_path, program_path, parameters_path, network_index=None, normalization=False):
            retry_count = 0
            backoff_time = 0.1
            attempts = 20
            result = os.system(f'{slurm_path} {program_path} {parameters_path}')

            while result != 0:
                retry_count += 1
                if retry_count > attempts:
                    if normalization:
                        print(
                            f"Failed to submit normalization job for network {network_index} after {attempts} attempts. Skipping.")
                    else:
                        print(f"Failed to submit job for network {network_index} after {attempts} attempts. Skipping.")
                    break

                if normalization:
                    print(
                        f"Retry #{retry_count} for normalization run network {network_index}, waiting {backoff_time} seconds...")
                else:
                    print(f"Retry #{retry_count} for network {network_index}, waiting {backoff_time} seconds...")

                time.sleep(backoff_time)
                backoff_time = min(backoff_time * 2, 1)
                result = os.system(f'{slurm_path} {program_path} {parameters_path}')

        # dir_path = os.path.dirname(os.path.realpath(__file__))
        # slurm_path = dir_path + '/slurm.serjob'
        # program_path = dir_path + '/cwesis.exe'
        # os.mkdir(foldername)
        # os.chdir(foldername)
        # data_path = os.getcwd() + '/'
        N, sims, start, k, x, lam, duartion, Num_inf, Alpha, number_of_networks, tau, eps_din, eps_dout, strength, prog, Beta_avg, error_graphs, correlation = parameters
        N, sims, start, k, x, lam, duartion, Num_inf, Alpha, number_of_networks, tau, \
        eps_din, eps_dout, strength, prog, Beta_avg, error_graphs, correlation = \
            int(N), int(sims), float(start), int(k), float(x), float(lam), float(duartion), int(Num_inf), float(
                Alpha), int(number_of_networks), float(tau), float(eps_din), float(eps_dout), \
            float(strength), prog, float(Beta_avg), bool(error_graphs), float(correlation)
        graph_degrees = np.array([G.degree(n) for n in G.nodes()])
        k_avg_graph,graph_std,graph_skewness = np.mean(graph_degrees),np.std(graph_degrees),skew(graph_degrees)
        second_moment,third_moment = np.mean((graph_degrees)**2),np.mean((graph_degrees)**3)
        eps_graph = graph_std / k_avg_graph
        largest_eigenvalue, largest_eigen_vector = eigsh(nx.adjacency_matrix(G).astype(float), k=1, which='LA',
                                                         return_eigenvectors=True)
        Beta = float(lam) / largest_eigenvalue[0]
        graph_correlation = nx.degree_assortativity_coefficient(G)
        Istar = (1 - 1 / lam) * N
        parameters = np.array(
            [N, sims, start, k_avg_graph, x, lam, Alpha, Beta, tau, Istar, strength, prog,
             dir_path, eps_graph, eps_graph, duartion, strength * Beta, graph_std, graph_skewness, third_moment, second_moment,graph_correlation])
        np.save('parameters_all.npy', parameters)
        infile = 'GNull.pickle'
        # correlation_graph = 2.0 if correlation != 0 else 0.0
        with open(infile, 'wb') as f:
            pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
        submit_job(N, prog, lam, eps_din, eps_dout, correlation, number_of_networks, k,
                   error_graphs, sims, tau, start, duartion, strength, relaxation_time, x,
                   Alpha, run_mc_simulation, normalization_run_flag, slurm_path, program_path,
                   correlation_heat_map,infile)


    if heatmap:
        # measurements = np.where(eps_din < 0.2, 1000000, 10000)
        measurements = np.where(correlation < 0.2, 1000000, 10000)
        sims = (measurements/number_of_networks).astype(int)
        loop_over = correlation
        duartion = np.linspace(0.01, 3.0, 2)
        if correlation_heat_map:
            G = rand_networks.configuration_model_undirected_graph_mulit_type(k, eps_din, N, prog, 0.0)
        for d in duartion:
            for i, j in zip(loop_over, sims):
                if correlation_heat_map==False:
                    submit_job(N, prog, lam, eps_din, eps_dout, i, number_of_networks, k,error_graphs, j, tau, start, d,
                               strength, relaxation_time, x,Alpha, run_mc_simulation, normalization_run_flag, slurm_path,
                               program_path)
                else:
                    Num_inf = int(x * N)
                    Beta_avg = Alpha * lam / k
                    parameters = np.array(
                        [N, j, start, k, x, lam, d, Num_inf, Alpha, number_of_networks, tau, eps_din,
                         eps_dout, strength, prog, Beta_avg, error_graphs, i])
                    submit_correlation_heatmap(G,correlation_heat_map,parameters)

    else:
        # measurements = 1000000
        # duartion = np.linspace(0.01,2.0,20)
        # sims = int(measurements/number_of_networks)
        loop_over = correlation
        measurements = np.where(loop_over < 1.0, 10000000, 1000000)
        sims = (measurements/number_of_networks).astype(int)
        for i,j in zip(loop_over,sims):
            submit_job(N, prog, lam, eps_din, eps_dout, i, number_of_networks, k,error_graphs, j, tau,
                       start, duartion, strength, relaxation_time, x,Alpha, run_mc_simulation, normalization_run_flag,
                       slurm_path, program_path)