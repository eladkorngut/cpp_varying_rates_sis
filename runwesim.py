import numpy as np
import os
import rand_networks
import csv
import pickle
import networkx as nx
from scipy.stats import skew
import argparse
import time
from scipy.sparse.linalg import eigsh

def export_parameters_to_csv(parameters,network_number):
    name_parameters = 'cparameters_{}.txt'.format(network_number)
    # N, sims, it, k, x, lam, jump, Alpha,Beta,number_of_networks, tau, mf_solution ,eps_din, eps_dout, new_trajcetory_bin, prog, Beta_avg,dir_path = parameters
    # cparameters=[N, sims, it, k, x, lam, jump, Alpha,Beta,number_of_networks, tau, mf_solution ,new_trajcetory_bin, prog, Beta_avg,dir_path]
    f =open(name_parameters,'+a')
    with f:
        writer = csv.writer(f)
        writer.writerow(parameters)
    f.close()

def export_network_to_csv(G,netname):
    # Open a CSV file for writing incoming neighbors

    # Check if the graph is directed
    is_directed = G.is_directed()

    with open('Adjin_{}.txt'.format(netname), 'w', newline='') as incoming_file:
        # Create a CSV writer
        incoming_writer = csv.writer(incoming_file)
        # Iterate over all nodes in the graph
        for node in np.sort(G):
            # Get the incoming neighbors of the current node
            if is_directed:
                incoming_neighbors = list(G.predecessors(node))
                # Get the degree of the current node
                degree = G.in_degree[node]
            else:
                incoming_neighbors = list(G.neighbors(node))  # All neighbors for undirected graph
                degree = G.degree[node]
            # Write a row to the CSV file for the current node
            joint = np.concatenate(([degree],incoming_neighbors),axis=0)
            incoming_writer.writerow(joint)
            # incoming_writer.writerow([degree])
            # for node in incoming_neighbors: incoming_writer.writerow([node])
    with open('Adjout_{}.txt'.format(netname), 'w', newline='') as outgoing_file:
        # Create a CSV writer
        outgoing_writer = csv.writer(outgoing_file)
        # Iterate over all nodes in the graph
        for node in np.sort(G):
            # Get the incoming neighbors of the current node
            if is_directed:
                outgoing_neighbors = list(G.successors(node))
                # Get the degree of the current node
                degree = G.out_degree[node]
            else:
                outgoing_neighbors = list(G.neighbors(node))  # All neighbors for undirected graph
                degree = G.degree[node]
            # Write a row to the CSV file for the current node
            joint = np.concatenate(([degree],outgoing_neighbors),axis=0)
            outgoing_writer.writerow(joint)

def job_to_cluster(foldername,parameters,Istar,normalization_run,runheatcorrelation,graphname):
    # This function submit jobs to the cluster with the following program keys:
    # bd: creates a bimodal directed networks and find its mean time to extinction

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

    dir_path = os.path.dirname(os.path.realpath(__file__))
    slurm_path = dir_path +'/slurm.serjob'
    program_path = dir_path +'/cwesis.exe'
    if runheatcorrelation:
        with open(graphname, 'rb') as f:
            G = pickle.load(f)
        graph_correlation = nx.degree_assortativity_coefficient(G)
    os.makedirs(foldername, exist_ok=True)
    # os.mkdir(foldername)
    os.chdir(foldername)
    data_path = os.getcwd() +'/'
    N, sims, start, k, x, lam, duartion, Num_inf, Alpha, number_of_networks, tau, eps_din, eps_dout, strength, prog, Beta_avg, error_graphs,correlation = parameters
    N, sims, start, k, x, lam, duartion, Num_inf, Alpha, number_of_networks, tau, \
    eps_din, eps_dout, strength, prog, Beta_avg, error_graphs,correlation=\
    int(N),int(sims),float(start),float(k),float(x),float(lam),float(duartion),int(Num_inf),float(Alpha),int(number_of_networks),float(tau),float(eps_din),float(eps_dout),\
    float(strength),prog,float(Beta_avg),bool(error_graphs),float(correlation)
    if 'prog'!='1d' and runheatcorrelation==False:
        G = rand_networks.configuration_model_undirected_graph_mulit_type(k,eps_din,N,prog,correlation)
        graph_degrees = np.array([G.degree(n) for n in G.nodes()])
        k_avg_graph,graph_std,graph_skewness = np.mean(graph_degrees),np.std(graph_degrees),skew(graph_degrees)
        second_moment,third_moment = np.mean((graph_degrees)**2),np.mean((graph_degrees)**3)
        eps_graph = graph_std / k_avg_graph
        largest_eigenvalue, largest_eigen_vector = eigsh(nx.adjacency_matrix(G).astype(float), k=1, which='LA',
                                                         return_eigenvectors=True)
        Beta = float(lam) / largest_eigenvalue[0]
        graph_correlation = nx.degree_assortativity_coefficient(G)
        parameters = np.array(
            [N, sims, start, k_avg_graph, x, lam, Alpha, Beta, tau, Istar, strength, prog,
             dir_path, eps_graph, eps_graph, duartion, strength * Beta, graph_std, graph_skewness, third_moment, second_moment,graph_correlation])
        np.save('parameters_all.npy', parameters)
    elif runheatcorrelation:
        G, graph_correlation = rand_networks.xulvi_brunet_sokolov_target_assortativity(G, correlation,graph_correlation, 0.05, 1000000)
        graph_degrees = np.array([G.degree(n) for n in G.nodes()])
        k_avg_graph,graph_std,graph_skewness = np.mean(graph_degrees),np.std(graph_degrees),skew(graph_degrees)
        second_moment,third_moment = np.mean((graph_degrees)**2),np.mean((graph_degrees)**3)
        eps_graph = graph_std / k_avg_graph
        largest_eigenvalue, largest_eigen_vector = eigsh(nx.adjacency_matrix(G).astype(float), k=1, which='LA',
                                                         return_eigenvectors=True)
        Beta = float(lam) / largest_eigenvalue[0]
        graph_correlation = nx.degree_assortativity_coefficient(G)
        parameters = np.array(
            [N, sims, start, k_avg_graph, x, lam, Alpha, Beta, tau, Istar, strength, prog,
             dir_path, eps_graph, eps_graph, duartion, strength * Beta, graph_std, graph_skewness, third_moment, second_moment,graph_correlation])
        np.save('parameters_all.npy', parameters)
    else:
        G = nx.random_regular_graph(k, N)
        eps_graph = 0
        k_avg_graph,graph_std,graph_skewness = k,0.0,0.0
        second_moment,third_moment = 0.0,0.0
        graph_correlation = 0.0
        Beta_graph = float(lam)/k_avg_graph
        Beta = Beta_graph
        parameters = np.array(
            [N, sims, start, k_avg_graph, x, lam, Alpha, Beta, tau, Istar, strength, prog,
             dir_path, eps_graph, eps_graph, duartion, strength * Beta, graph_std, graph_skewness, third_moment, second_moment,graph_correlation])
        np.save('parameters_all.npy', parameters)

    for i in range(int(number_of_networks)):
        if error_graphs==False:
            if prog!='1d' and runheatcorrelation==False:
                G = rand_networks.configuration_model_undirected_graph_mulit_type(float(k),float(eps_din),int(N),prog)
                graph_degrees = np.array([G.degree(n) for n in G.nodes()])
                k_avg_graph, graph_std, graph_skewness = np.mean(graph_degrees), np.std(graph_degrees), skew(graph_degrees)
                second_moment,third_moment = np.mean((graph_degrees)**2),np.mean((graph_degrees)**3)
                eps_graph = graph_std / k_avg_graph
                # third_moment = graph_skewness * (graph_std ** 3)
                largest_eigenvalue, largest_eigen_vector = eigsh(nx.adjacency_matrix(G).astype(float), k=1, which='LA',
                                                                 return_eigenvectors=True)
                Beta = float(lam) / largest_eigenvalue[0]
                infile = 'GNull_{}.pickle'.format(i)
                with open(infile, 'wb') as f:
                    pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
                nx.write_gpickle(G, infile)
            elif runheatcorrelation:
                with open(graphname, 'rb') as f:
                    G = pickle.load(f)
                G, correlation_graph = rand_networks.xulvi_brunet_sokolov_target_assortativity(G, correlation_factor,
                                    correlation_graph, 0.05,1000000)
                graph_degrees = np.array([G.degree(n) for n in G.nodes()])
                k_avg_graph, graph_std, graph_skewness = np.mean(graph_degrees), np.std(graph_degrees), skew(
                    graph_degrees)
                second_moment, third_moment = np.mean((graph_degrees) ** 2), np.mean((graph_degrees) ** 3)
                eps_graph = graph_std / k_avg_graph
                largest_eigenvalue, largest_eigen_vector = eigsh(nx.adjacency_matrix(G).astype(float), k=1, which='LA',
                                                                 return_eigenvectors=True)
                Beta = float(lam) / largest_eigenvalue[0]
                infile = 'GNull_{}.pickle'.format(i)
                with open(infile, 'wb') as f:
                    pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
                graph_correlation = nx.degree_assortativity_coefficient(G)
                parameters = np.array(
                    [N, sims, start, k_avg_graph, x, lam, Alpha, Beta, tau, Istar, strength, prog,
                     dir_path, eps_graph, eps_graph, duartion, strength * Beta, graph_std, graph_skewness, third_moment,
                     second_moment, graph_correlation])
                np.save('parameters_all.npy', parameters)
            else:
                G = nx.random_regular_graph(k,N)
                k_avg_graph, graph_std, graph_skewness = k, 0.0, 0.0
                second_moment,third_moment = 0.0,0.0
                eps_graph = graph_std / k_avg_graph
        graph_correlation = nx.degree_assortativity_coefficient(G)
        parameters = np.array([N,sims,start,k_avg_graph,x,lam,Alpha,Beta,i,tau,Istar,strength,prog,dir_path,eps_graph,
                               eps_graph,duartion,strength*Beta,graph_std,graph_skewness,third_moment,second_moment,graph_correlation])
        np.save('parameters_{}.npy'.format(i), parameters)
        # infile = 'GNull_{}.pickle'.format(i)
        # with open(infile,'wb') as f:
        #     pickle.dump(G,f,pickle.HIGHEST_PROTOCOL)
        # nx.write_gpickle(G, infile)
        export_network_to_csv(G, i)
        export_parameters_to_csv(parameters,i)
        path_adj_in = data_path + 'Adjin_{}.txt'.format(i)
        path_adj_out = data_path + 'Adjout_{}.txt'.format(i)
        path_parameters = data_path + 'cparameters_{}.txt'.format(i)
        parameters_path = '{} {} {}'.format(path_adj_in,path_adj_out,path_parameters)
        submit_with_retries(slurm_path, program_path, parameters_path, network_index=i)

        if normalization_run:
            # Go back to the parent directory
            os.chdir('..')
            # Create a new folder for the normalization run; for example, append '_normalization' to the original folder name.
            norm_folder = foldername + '_normalization'
            if not os.path.exists(norm_folder):
                os.mkdir(norm_folder)
            # Change to the newly created normalization folder.
            os.chdir(norm_folder)
            parameters_normalization = np.array([
                N, sims, start, k_avg_graph, x, lam, Alpha, Beta, i, tau, Istar,
                1.0, prog, dir_path, eps_graph, eps_graph, duartion, Beta,
                graph_std, graph_skewness, third_moment, second_moment,graph_correlation])
            parameters_normalization_all = np.array([
                N, sims, start, k_avg_graph, x, lam, Alpha, Beta, tau, Istar,
                1.0, prog, dir_path, eps_graph, eps_graph, duartion, Beta,
                graph_std, graph_skewness, third_moment, second_moment,graph_correlation])
            np.save('parameters_all.npy', parameters_normalization_all)
            np.save('parameters_{}.npy'.format(i), parameters_normalization)
            export_network_to_csv(G, i)
            export_parameters_to_csv(parameters_normalization, i)
            data_path_norm = os.getcwd() + '/'
            np.save('parameters_{}.npy'.format(i), parameters_normalization)
            path_adj_in_norm = data_path_norm + 'Adjin_{}.txt'.format(i)
            path_adj_out_norm = data_path_norm + 'Adjout_{}.txt'.format(i)
            path_parameters_norm = data_path_norm + 'cparameters_{}.txt'.format(i)
            parameters_path_norm = '{} {} {}'.format(path_adj_in_norm,path_adj_out_norm,path_parameters_norm)
            submit_with_retries(slurm_path, program_path, parameters_path_norm, network_index=i, normalization=True)
            os.chdir(data_path)
        # os.system('{} {} {} {}'.format(program_path,path_adj_in,path_adj_out,path_parameters))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process network and WE method parameters.")

    # Parameters for the network
    parser.add_argument('--N', type=int, help='Number of nodes')
    parser.add_argument('--prog', type=str, help='Program')
    parser.add_argument('--lam', type=float, help='The reproduction number')
    parser.add_argument('--eps_din', type=float, help='The normalized std (second moment divided by the first) of the in-degree distribution')
    parser.add_argument('--eps_dout', type=float, help='The normalized std (second moment divided by the first) of the out-degree distribution')
    parser.add_argument('--correlation', type=float, help='Correlation parameter')
    parser.add_argument('--number_of_networks', type=int, help='Number of networks')
    parser.add_argument('--k', type=int, help='Average number of neighbors for each node')
    parser.add_argument('--error_graphs', action='store_true', help='Flag for error graphs')

    # Parameters for the WE method
    parser.add_argument('--sims', type=int, help='Number of simulations at each bin')
    parser.add_argument('--tau', type=int, help='Tau parameter')
    parser.add_argument('--start', type=int, help='The time for which the catastrophe starts')
    parser.add_argument('--duartion', type=float, help='The duration, how long does the catstrophe lasts')
    parser.add_argument('--strength', type=float, help='By how much the decline in Beta, how strong is the catastrophe')

    # Parameters that don't get changed
    parser.add_argument('--relaxation_time', type=int, help='Relaxation time')
    parser.add_argument('--x', type=float, help='Initial infection percentage')
    parser.add_argument('--Alpha', type=float, help='Recovery rate')
    parser.add_argument('--run_mc_simulation', action='store_true', help='Flag to run MC simulation')
    parser.add_argument('--short_path', action='store_true', help='Flag to measure mean shortest path')
    parser.add_argument('--normalization_run', action='store_true', help='Flag for running parameter normalization')
    parser.add_argument('--runheatcorrelation', action='store_true', help='Flag for xuli-soklov graph')
    parser.add_argument('--graphname', type=str, help='graph name')


    args = parser.parse_args()
    normalization_run = args.normalization_run

    # Default parameters
    N = 1000 if args.N is None else args.N
    prog = 'gam' if args.prog is None else args.prog
    graphname = 'GNull.pickle' if args.graphname is None else args.graphname
    lam = 1.2 if args.lam is None else args.lam
    eps_din = 0.5 if args.eps_din is None else args.eps_din
    eps_dout = 0.5 if args.eps_dout is None else args.eps_dout
    correlation = -0.1 if args.correlation is None else args.correlation
    number_of_networks = 10 if args.number_of_networks is None else args.number_of_networks
    k = 50 if args.k is None else args.k
    error_graphs = args.error_graphs
    normalization_run = args.normalization_run
    # runheatcorrelation = args.runheatcorrelation
    runheatcorrelation = True
    # normalization_run = True


    sims = 1000 if args.sims is None else args.sims
    tau = 150 if args.tau is None else args.tau
    start = 80 if args.start is None else args.start
    duartion = 1.0 if args.duartion is None else args.duartion
    strength = 0.0 if args.strength is None else args.strength

    relaxation_time = 20 if args.relaxation_time is None else args.relaxation_time
    x = 0.2 if args.x is None else args.x
    Num_inf = int(x * N)
    Alpha = 1.0 if args.Alpha is None else args.Alpha
    Beta_avg = Alpha * lam / k
    # run_mc_simulationtion = True


    # parameters = np.array([N, sims, start, k, x, lam, duartion, Num_inf, Alpha, number_of_networks, tau, eps_din,
    #                        eps_dout, strength, prog, Beta_avg, error_graphs,correlation])

    parameters = [
        int(N),                   # --N
        int(sims),                # --sims
        int(start),               # --start
        int(k),                   # --k
        float(x),                 # --x
        float(lam),               # --lam
        float(duartion),          # --duartion
        int(Num_inf),             # seed
        float(Alpha),             # --Alpha
        int(number_of_networks),  # --number_of_networks
        float(tau),               # --tau
        float(eps_din),           # --eps_din
        float(eps_dout),          # --eps_dout
        float(strength),          # --strength
        prog,                     # --prog  (string)
        float(Beta_avg),          # computed Beta
        int(error_graphs),        # 0 or 1
        float(correlation)        # --correlation
    ]



    foldername = 'prog_{}_N{}_k_{}_R_{}_tau_{}_start_{}_duartion_{}_strength_{}_sims_{}_net_{}_epsin_{}_epsout_{}_correlation_{}_err_{}'.format(
        prog, N, k, lam, tau, start, duartion, strength, sims, number_of_networks, eps_din, eps_dout, correlation,error_graphs)
    Istar = (1 - 1/lam) * N
    job_to_cluster(foldername, parameters, Istar,normalization_run,runheatcorrelation,graphname)
    # act_as_main(foldername, parameters, Istar, prog)