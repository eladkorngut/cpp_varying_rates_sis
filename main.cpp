#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <random>
#include <sstream>
#include <fstream>
#include <string>
#include <tuple>
#include <numeric>
#include <experimental/filesystem>
#include <deque>
#include <tr1/tuple>
#include <chrono>

// Declaration of functions
void remove_infected_node(int node,std::vector<int>& infected_node,std::vector<int> &positions);
void add_susceptible_node(int node,int m,std::vector<std::vector<int>>& susceptible_node,std::vector<int>& positions);
void decrement_susc_neighbs(int k_in,int k_out,std::vector<int> &neighbs_in,std::vector<int> &neighbs_out,
                            std::vector<std::vector<int>> &net_susceptible_nodes,std::vector<int> &infected_neighbors_in,
                            std::vector<int> &infected_neighbors_out,std::vector<int> &s_m,std::vector<int> &positions
        ,std::vector<int> &sigma,std::vector<int> &inf_neigh);
void increment_susc_neighbs(int k_in,int k_out,std::vector<int> &neighbs_in,std::vector<int> &neighbs_out,
                            std::vector<std::vector<int>> &net_suceptible_nodes,std::vector<int> &infected_neighbors_in,
                            std::vector<int> &infected_neighbors_out,std::vector<int> &s_m,std::vector<int> &positions,
                            std::vector<int> &sigma);
void remove_susceptible_node(int node,int k,std::vector<std::vector<int>>& susceptible_nodes,std::vector<int>& positions);
void add_infected_node(int node,std::vector<int>& infected_node,std::vector<int>& positions);
void read_in_neighborslist(int N,std::string& filename,std::vector<std::vector<int>>& Adjlist,std::vector<int>& degree,std::fstream& Adjfile);
void read_parameters(std::string& filename,std::fstream& parameters_file,int N,int sims,int it,int k,double x,double lam,double Alpha,double Beta,int network_number,
                     double mf_solution,int new_trajectory_bin,int k_max);

std::pair<double,int> gillespie1d(double tau,std::mt19937 &gen,std::uniform_real_distribution<double>& uniform_dist,
                   double cat_start,double cat_duration,int N, int seed,double Alpha, double beta, double beta_cat,int k){
    double time = 0.0,rates,r1,r2;
    int Num_inf = seed;
    while (Num_inf>0 && time<tau){
        r1 = uniform_dist(gen);
        r2 = -std::log(uniform_dist(gen));
        if (time > cat_start and time <= cat_start + cat_duration)
            rates = Num_inf * Alpha + (beta_cat * k * (N - Num_inf) * Num_inf)/N;
        else
            rates = Num_inf * Alpha + (beta * k * (N - Num_inf) * Num_inf)/N;
        if (Alpha*Num_inf/rates<r1)
            Num_inf++;
        else
            Num_inf--;
        time+= r2 / rates;
    }
//    if (Num_inf>0.0)
//        return 0.0;
    std::pair<double,int> network_data(time,Num_inf);
    return network_data;
}


struct network_topology{
public:
    double Alpha,Beta,beta_cat;
    int k_max;
    std::vector<int> degree_in,degree_out;
    std::vector<std::vector<int>> Adjlist_in,Adjlist_out;

    network_topology(double Alpha,double Beta,int k_max,std::vector<int> &degree_in,std::vector<int> &degree_out,
                     std::vector<std::vector<int>> &Adjlist_in,std::vector<std::vector<int>> &Adjlist_out,double beta_cat): Alpha(Alpha),
                                                                                                            Beta(Beta),k_max(k_max),degree_in(degree_in),degree_out(degree_out),Adjlist_in(Adjlist_in),Adjlist_out(Adjlist_out),beta_cat(beta_cat){}
    void decrement_susc_neighbs(int k_in,int k_out,std::vector<int> &neighbs_in,std::vector<int> &neighbs_out,
                                std::vector<std::vector<int>> &net_susceptible_nodes,std::vector<int> &infected_neighbors_in,
                                std::vector<int> &infected_neighbors_out,std::vector<int> &s_m,std::vector<int> &positions
            ,std::vector<int> &sigma,std::vector<int> &inf_neigh){
        int neighb,m2;
        for (int j=0;j<k_out;j++){
            // Change neighbours from s_m class to s_m{m-1} class. need to be for the ingoing neighbours
            neighb = neighbs_out[j];
            infected_neighbors_in[neighbs_out[j]]--;
            if (sigma[neighb]==0){
                m2 = infected_neighbors_in[neighb];
                remove_susceptible_node(neighb,m2+1,net_susceptible_nodes,positions);
                s_m[m2+1]--;
                add_susceptible_node(neighb,m2,net_susceptible_nodes,positions);
                s_m[m2]++;
            }

        }
        for (int j=0;j<k_in;j++){
            infected_neighbors_out[neighbs_in[j]]--; // Each neighbor of recovered node has one less infected neighbor
        }
    }

};

struct networks_dynamics{
public:
    std::list<int> num_inf;
    std::list<double> weights;
    std::list<double> avec_sum;
    std::list<double> t;
    std::list<std::vector<int>> infected_node,infected_neighbors_in,infected_neighbors_out,sigma,s_m,positions;
    std::list<std::vector<std::vector<int>>> susceptible_nodes;
    std::list<int> SI;


    // Constructor to initialize the member variables
    networks_dynamics(std::list<int> &num_inf,std::list<double> &weights,std::list<double> &avec_sum,
                      std::list<double> &t,std::list<std::vector<int>> &infected_node,std::list<std::vector<int>> &infected_neighbors_in,
                      std::list<std::vector<int>> &infected_neighbors_out,std::list<std::vector<int>> &sigma,std::list<std::vector<int>> &s_m,
                      std::list<std::vector<int>> &positions,std::list<std::vector<std::vector<int>>> &susceptible_nodes,std::list<int> &SI):
            num_inf(num_inf),weights(weights),avec_sum(avec_sum),infected_node(infected_node),infected_neighbors_in(infected_neighbors_in),infected_neighbors_out(infected_neighbors_out),
            sigma(sigma),positions(positions),susceptible_nodes(susceptible_nodes),SI(SI),s_m(s_m),t(t){}

};


struct simulation_data{
public:
    std::list<int>::iterator i_num_inf,i_SI;
    std::list<double>::iterator i_weights,i_avec_sum,i_t;
    std::list<std::vector<int>>::iterator i_infected_node,i_infected_neighbors_in,i_infected_neighbors_out,i_sigma,i_s_m,i_positions;
    std::list<std::vector<std::vector<int>>>::iterator i_susceptible_nodes;

    // Constructor to initialize the member variables
    simulation_data(
            // The types are: num_inf,weights, avec_sum,t,infected_node,infected_neighbors_in,infected_neighbors_out,sigma,s_m,positions,susceptible_nodes,SI
            std::list<int>::iterator i_num_inf,
            std::list<double>::iterator i_weights,
            std::list<double>::iterator i_avec_sum,
            std::list<double>::iterator i_t,
            std::list<std::vector<int>>::iterator i_infected_node,
            std::list<std::vector<int>>::iterator i_infected_neighbors_in,
            std::list<std::vector<int>>::iterator i_infected_neighbors_out,
            std::list<std::vector<int>>::iterator i_sigma,
            std::list<std::vector<int>>::iterator i_s_m,
            std::list<std::vector<int>>::iterator i_positions,
            std::list<std::vector<std::vector<int>>>::iterator i_susceptible_nodes,
            std::list<int>::iterator i_SI
    ) : i_num_inf(i_num_inf), i_weights(i_weights), i_avec_sum(i_avec_sum), i_t(i_t), i_infected_node(i_infected_node), i_infected_neighbors_in(i_infected_neighbors_in)
            , i_infected_neighbors_out(i_infected_neighbors_out), i_sigma(i_sigma), i_s_m(i_s_m), i_positions(i_positions), i_susceptible_nodes(i_susceptible_nodes), i_SI(i_SI) {}
    simulation_data(
            // The types are: num_inf,weights, avec_sum,t,infected_node,infected_neighbors_in,infected_neighbors_out,sigma,s_m,positions,susceptible_nodes,SI
            std::list<int> &num_inf,
            std::list<double> &weights,
            std::list<double> &avec_sum,
            std::list<double> &t,
            std::list<std::vector<int>> &infected_node,
            std::list<std::vector<int>> &infected_neighbors_in,
            std::list<std::vector<int>> &infected_neighbors_out,
            std::list<std::vector<int>> &sigma,
            std::list<std::vector<int>> &s_m,
            std::list<std::vector<int>> &positions,
            std::list<std::vector<std::vector<int>>> &susceptible_nodes,
            std::list<int> &SI
    ) : i_num_inf(num_inf.begin()), i_weights(weights.begin()), i_avec_sum(avec_sum.begin()), i_t(t.begin()), i_infected_node(infected_node.begin()), i_infected_neighbors_in(infected_neighbors_in.begin())
            , i_infected_neighbors_out(infected_neighbors_out.begin()), i_sigma(sigma.begin()), i_s_m(s_m.begin()), i_positions(positions.begin()), i_susceptible_nodes(susceptible_nodes.begin()), i_SI(SI.begin()) {}

    simulation_data (networks_dynamics& net_d):i_num_inf(net_d.num_inf.begin()),i_weights(net_d.weights.begin()),i_avec_sum(net_d.avec_sum.begin()),i_t(net_d.t.begin()),i_infected_node(net_d.infected_node.begin()),
                                               i_infected_neighbors_in(net_d.infected_neighbors_in.begin()),i_infected_neighbors_out(net_d.infected_neighbors_out.begin()),i_sigma(net_d.sigma.begin()),i_s_m(net_d.s_m.begin()),i_positions(net_d.positions.begin()),
                                               i_susceptible_nodes(net_d.susceptible_nodes.begin()),i_SI(net_d.SI.begin()){}


    simulation_data& operator++(){
        ++i_num_inf;
        ++i_weights;
        ++i_avec_sum;
        ++i_t;
        ++i_infected_node;
        ++i_infected_neighbors_in;
        ++i_infected_neighbors_out;
        ++i_sigma;
        ++i_s_m;
        ++i_positions;
        ++i_susceptible_nodes;
        ++i_SI;
        return *this;
    }
    bool end(networks_dynamics &net_d){
        if (net_d.num_inf.end()==i_num_inf){
            return true;
        }
        return false;
    }

    simulation_data erase_simulation(networks_dynamics &net_d){
        this->i_num_inf = net_d.num_inf.erase(i_num_inf);
        this->i_weights = net_d.weights.erase(i_weights);
        this->i_avec_sum = net_d.avec_sum.erase(i_avec_sum);
        this->i_t = net_d.t.erase(i_t);
        this->i_infected_node = net_d.infected_node.erase(i_infected_node);
        this->i_infected_neighbors_in = net_d.infected_neighbors_in.erase(i_infected_neighbors_in);
        this->i_infected_neighbors_out = net_d.infected_neighbors_out.erase(i_infected_neighbors_out);
        this->i_sigma = net_d.sigma.erase(i_sigma);
        this->i_s_m = net_d.s_m.erase(i_s_m);
        this->i_positions = net_d.positions.erase(i_positions);
        this->i_susceptible_nodes = net_d.susceptible_nodes.erase(i_susceptible_nodes);
        this->i_SI = net_d.SI.erase(i_SI);
        return *this;
    }

    std::tuple<double,int,simulation_data> gillespie(double tau,std::mt19937 &gen,std::uniform_real_distribution<double>& uniform_dist,
                                                 networks_dynamics &net_d,network_topology &net_t,double cat_start,double cat_duration){
        double sim_time=0,r1,r2,s_m_sum=0;
        int m_in,m_out,node,i_num_inf_previous;
        while(*i_num_inf>0 && sim_time<tau){
            r1 = uniform_dist(gen);
            r2 = -std::log(uniform_dist(gen));
            if (sim_time>cat_start && sim_time<=cat_start+cat_duration)
                *i_avec_sum = (*i_num_inf)*net_t.Alpha+(*i_SI)*net_t.beta_cat;
            else
                *i_avec_sum = (*i_num_inf)*net_t.Alpha+(*i_SI)*net_t.Beta;
//            *i_avec_sum = (*i_num_inf)*net_t.Alpha+(*i_SI)*net_t.Beta;
            sim_time+=r2/(*i_avec_sum);
            //Pick a node, change its state and update the transition rates
            if (r1<(*i_num_inf)*net_t.Alpha/(*i_avec_sum)) {
                // Recover an infected node
                std::uniform_int_distribution<int> node_dist(0,*i_num_inf-1);
                node = (*i_infected_node)[node_dist(gen)]; //randomly choose infected node
                m_in = (*i_infected_neighbors_in)[node]; //number of infected neighbors that infecting the node
                m_out = (*i_infected_neighbors_out)[node]; //number of infected neighbors that the node infect
                (*i_sigma)[node]=0;
                (*i_num_inf)--;
                //Remove node from infected list and add to susceptible list
                remove_infected_node(node,*i_infected_node,*i_positions);
                (*i_SI) = (*i_SI) + m_in + m_out - net_t.degree_out[node];
                add_susceptible_node(node,m_in,*i_susceptible_nodes,*i_positions);
                // Increment number of susceptible nodes with m infected nodes
                (*i_s_m)[m_in]++;
                decrement_susc_neighbs(net_t.degree_in[node],net_t.degree_out[node],net_t.Adjlist_in[node],net_t.Adjlist_out[node],*i_susceptible_nodes,
                                       *i_infected_neighbors_in,*i_infected_neighbors_out,*i_s_m,*i_positions,*i_sigma,*i_infected_node);
            }
            else{
                //Infect a susceptible node
                // Pick an m class with probability proportional to the total number of SI links in that m class i.e. s_m[m]*m
                r1= uniform_dist(gen)*(*i_SI);
                s_m_sum = 0;
                for (int m1=1;m1<=net_t.k_max;m1++){
                    s_m_sum = s_m_sum + (*i_s_m)[m1]*m1;
                    if (r1<=s_m_sum){
                        // choose a node with m infected neighbours at random
                        std::uniform_int_distribution<int> node_dist(0,(*i_s_m)[m1]-1);
                        node = (*i_susceptible_nodes)[m1][node_dist(gen)];
                        break;
                    }
                }
                m_in = (*i_infected_neighbors_in)[node];
                m_out = (*i_infected_neighbors_out)[node];
                (*i_sigma)[node] = 1;
                *i_num_inf = (*i_num_inf)+1;
                // Remove the node from susceptible list and add to infected list
                remove_susceptible_node(node,m_in,*i_susceptible_nodes,*i_positions);
                add_infected_node(node,*i_infected_node,*i_positions);
                //Decrement number of susceptible nodes
                (*i_s_m)[m_in]--;
                // change the neighbours to other lists and adjust s_m
                (*i_SI) = (*i_SI) + net_t.degree_out[node] - m_in - m_out;
                increment_susc_neighbs(net_t.degree_in[node],net_t.degree_out[node],net_t.Adjlist_in[node],net_t.Adjlist_out[node],
                                       *i_susceptible_nodes,*i_infected_neighbors_in,*i_infected_neighbors_out,*i_s_m,*i_positions,*i_sigma);
            }
        }
        *i_t += sim_time;
        if (*i_num_inf<=0){
            i_num_inf_previous = *i_num_inf;
            return std::make_tuple(sim_time,i_num_inf_previous,erase_simulation(net_d));
        }
        return std::make_tuple(sim_time,*i_num_inf,*this);
    }
};

//        std::pair<double,simulation_data > gillespie(double tau,std::mt19937 &gen,std::uniform_real_distribution<double>& uniform_dist,
//                                                 networks_dynamics &net_d,network_topology &net_t,double cat_start,double cat_duration){
//        double sim_time=0,r1,r2,s_m_sum=0;
//        int m_in,m_out,node;
//        while(*i_num_inf>0 && sim_time<tau){
//            r1 = uniform_dist(gen);
//            r2 = -std::log(uniform_dist(gen));
//            if (sim_time>cat_start && sim_time<=cat_start+cat_duration)
//                *i_avec_sum = (*i_num_inf)*net_t.Alpha+(*i_SI)*net_t.beta_cat;
//            else
//                *i_avec_sum = (*i_num_inf)*net_t.Alpha+(*i_SI)*net_t.Beta;
////            *i_avec_sum = (*i_num_inf)*net_t.Alpha+(*i_SI)*net_t.Beta;
//            sim_time+=r2/(*i_avec_sum);
//            //Pick a node, change its state and update the transition rates
//            if (r1<(*i_num_inf)*net_t.Alpha/(*i_avec_sum)) {
//                // Recover an infected node
//                std::uniform_int_distribution<int> node_dist(0,*i_num_inf-1);
//                node = (*i_infected_node)[node_dist(gen)]; //randomly choose infected node
//                m_in = (*i_infected_neighbors_in)[node]; //number of infected neighbors that infecting the node
//                m_out = (*i_infected_neighbors_out)[node]; //number of infected neighbors that the node infect
//                (*i_sigma)[node]=0;
//                (*i_num_inf)--;
//                //Remove node from infected list and add to susceptible list
//                remove_infected_node(node,*i_infected_node,*i_positions);
//                (*i_SI) = (*i_SI) + m_in + m_out - net_t.degree_out[node];
//                add_susceptible_node(node,m_in,*i_susceptible_nodes,*i_positions);
//                // Increment number of susceptible nodes with m infected nodes
//                (*i_s_m)[m_in]++;
//                decrement_susc_neighbs(net_t.degree_in[node],net_t.degree_out[node],net_t.Adjlist_in[node],net_t.Adjlist_out[node],*i_susceptible_nodes,
//                                       *i_infected_neighbors_in,*i_infected_neighbors_out,*i_s_m,*i_positions,*i_sigma,*i_infected_node);
//            }
//            else{
//                //Infect a susceptible node
//                // Pick an m class with probability proportional to the total number of SI links in that m class i.e. s_m[m]*m
//                r1= uniform_dist(gen)*(*i_SI);
//                s_m_sum = 0;
//                for (int m1=1;m1<=net_t.k_max;m1++){
//                    s_m_sum = s_m_sum + (*i_s_m)[m1]*m1;
//                    if (r1<=s_m_sum){
//                        // choose a node with m infected neighbours at random
//                        std::uniform_int_distribution<int> node_dist(0,(*i_s_m)[m1]-1);
//                        node = (*i_susceptible_nodes)[m1][node_dist(gen)];
//                        break;
//                    }
//                }
//                m_in = (*i_infected_neighbors_in)[node];
//                m_out = (*i_infected_neighbors_out)[node];
//                (*i_sigma)[node] = 1;
//                *i_num_inf = (*i_num_inf)+1;
//                // Remove the node from susceptible list and add to infected list
//                remove_susceptible_node(node,m_in,*i_susceptible_nodes,*i_positions);
//                add_infected_node(node,*i_infected_node,*i_positions);
//                //Decrement number of susceptible nodes
//                (*i_s_m)[m_in]--;
//                // change the neighbours to other lists and adjust s_m
//                (*i_SI) = (*i_SI) + net_t.degree_out[node] - m_in - m_out;
//                increment_susc_neighbs(net_t.degree_in[node],net_t.degree_out[node],net_t.Adjlist_in[node],net_t.Adjlist_out[node],
//                                       *i_susceptible_nodes,*i_infected_neighbors_in,*i_infected_neighbors_out,*i_s_m,*i_positions,*i_sigma);
//            }
//        }
//        *i_t += sim_time;
//        if (*i_num_inf<=0){
//            double extinction_time=sim_time;
//            return std::make_pair(extinction_time,erase_simulation(net_d));
//        }
//        return std::make_pair(0.0,*this);
//    }
//};

// Start of auxiliary functions

int read_single_int_parameter(std::stringstream& ss,std::string& word){
    getline(ss, word, ',');
    ss.clear();
    return std::stoi(word);

}

double read_single_double_parameter(std::stringstream& ss,std::string& word) {
    getline(ss, word, ',');
    ss.clear();
    return std::stof(word);
}

std::string read_single_string_parameter(std::stringstream& ss,std::string& word) {
    getline(ss, word, ',');
    ss.clear();
    return word;
}

std::tuple<int,int,double,int,double,double,double,double,int,double,double,double,std::string & ,std::string &,double,double,double,double > read_parameters(std::string& filename,std::fstream& parameters_file) {
    // Read network parameters from python files
    std::string line;
    std::stringstream ss;
    std::string word;
    int N,sims,k,jump,network_number;
    double x,lam,Alpha,Beta,tau,mf_solution,duration,beta_cat,eps_din,eps_dout,strength,start;
    std::string prog,dir_path;
    parameters_file.open(filename);
    while (getline(parameters_file, line)) {
        ss.clear();
        ss.str(line);
        N = read_single_int_parameter(ss,word);
        sims = read_single_int_parameter(ss,word);
        start = read_single_int_parameter(ss,word);
        k = read_single_int_parameter(ss,word) ;
        x = read_single_double_parameter(ss,word);
        lam = read_single_double_parameter(ss,word);
//        jump = read_single_int_parameter(ss,word);
        Alpha = read_single_double_parameter(ss,word);
        Beta = read_single_double_parameter(ss,word);
        network_number = read_single_int_parameter(ss,word);
        tau = read_single_double_parameter(ss,word);
        mf_solution = read_single_double_parameter(ss,word);
        strength = read_single_int_parameter(ss,word);
        prog=read_single_string_parameter(ss,word);
        dir_path=read_single_string_parameter(ss,word);
        eps_din = read_single_double_parameter(ss,word);
        eps_dout = read_single_double_parameter(ss,word);
        duration = read_single_double_parameter(ss,word);
        beta_cat = read_single_double_parameter(ss,word);
    }
    return std::tuple<int,int,double,int,double,double,double,double,int,double,double,double,std::string & ,std::string&,double,double,double,double>
                    (N,sims,start,k,x,lam,Alpha,Beta,network_number,tau,mf_solution,strength,prog,dir_path,eps_din,eps_dout,duration,beta_cat);
}

void remove_infected_node(int node,std::vector<int> &infected_node,std::vector<int> &positions){
    // Remove the node from the list by swapping the node with the last node in the list and removing the last node in the list. This will not affect the positioning of the others
    infected_node[positions[node]] = infected_node.back();
    positions[infected_node.back()] = positions[node];
    infected_node.pop_back();
}

void add_susceptible_node(int node,int m,std::vector<std::vector<int>>& susceptible_node,std::vector<int>& positions){
    positions[node] = susceptible_node[m].size();
    susceptible_node[m].push_back(node);
}


void decrement_susc_neighbs(int k_in,int k_out,std::vector<int> &neighbs_in,std::vector<int> &neighbs_out,
                            std::vector<std::vector<int>> &net_susceptible_nodes,std::vector<int> &infected_neighbors_in,
                            std::vector<int> &infected_neighbors_out,std::vector<int> &s_m,std::vector<int> &positions
        ,std::vector<int> &sigma,std::vector<int> &inf_neigh){
    int neighb,m2;
    for (int j=0;j<k_out;j++){
        // Change neighbours from s_m class to s_m{m-1} class. need to be for the ingoing neighbours
        neighb = neighbs_out[j];
        infected_neighbors_in[neighbs_out[j]]--;
        if (sigma[neighb]==0){
            m2 = infected_neighbors_in[neighb];
            remove_susceptible_node(neighb,m2+1,net_susceptible_nodes,positions);
            s_m[m2+1]--;
            add_susceptible_node(neighb,m2,net_susceptible_nodes,positions);
            s_m[m2]++;
        }

    }
    for (int j=0;j<k_in;j++){
        infected_neighbors_out[neighbs_in[j]]--; // Each neighbor of recovered node has one less infected neighbor
    }
}

void remove_susceptible_node(int node,int m, std::vector<std::vector<int>> &net_susceptible_nodes,std::vector<int> &positions){

    // Remove the node from the list by swapping the node with the last node in the list. This will not affect the
    // positioning of the other nodes once the position of the swapped last node is accounted for
    net_susceptible_nodes[m][positions[node]] = net_susceptible_nodes[m].back();  // Move last node in list to nodes position
    positions[net_susceptible_nodes[m].back()] = positions[node]; // account for this change in the position vector
    net_susceptible_nodes[m].pop_back(); // remove the last node from the list
}

void increment_susc_neighbs(int k_in,int k_out,std::vector<int> &neighbs_in,std::vector<int> &neighbs_out,
                            std::vector<std::vector<int>> &net_suceptible_nodes,std::vector<int> &infected_neighbors_in,
                            std::vector<int> &infected_neighbors_out,std::vector<int> &s_m,std::vector<int> &positions,
                            std::vector<int> &sigma){
    int neighb,m2;
    for (int j=0;j<k_in;j++){
        infected_neighbors_out[neighbs_in[j]]++;
    }
    for (int j=0;j<k_out;j++){
        neighb = neighbs_out[j];
        infected_neighbors_in[neighb]++;
        if (sigma[neighb]==0){
            m2 = infected_neighbors_in[neighb];
            remove_susceptible_node(neighb,m2-1,net_suceptible_nodes,positions);
            s_m[m2-1]--;
            add_susceptible_node(neighb,m2,net_suceptible_nodes,positions);
            s_m[m2]++;
        }
    }
}


void add_infected_node(int node,std::vector<int> &infected_nodes,std::vector<int> &positions){
    positions[node]= infected_nodes.size();
    infected_nodes.push_back(node);
}


void read_in_neighborslist(int N,std::string& filename,std::vector<std::vector<int>> &Adjlist,std::vector<int> &degrees,std::fstream &Adjfile,int &k_max) {
    // Read network data from python files

    Adjfile.open(filename);
    std::string line;
    std::stringstream ss;
    int count_degree = 0;
    std::string word;
    if (Adjfile.is_open()) {
        while (std::getline(Adjfile, line)) {
            ss.clear();
            ss.str(line);
            std::getline(ss, word, ',');
            degrees[count_degree] = std::stoi(word);
            if (degrees[count_degree]>k_max)
                k_max = degrees[count_degree];
            while (getline(ss, word, ','))
                Adjlist[count_degree].push_back(std::stoi(word));
            count_degree++;
        }
    }
    Adjfile.close();
}

void increment_neighbours(int k,std::vector<int> &neighbs,std::vector<int> &infected_neighbs_in){
    for (int i=0;i<k;i++)
        infected_neighbs_in[neighbs[i]]++;
}

void intalized_network_random_infection(int N,int inital_infecteons,std::list< std::vector<int>>& sigma,std::list<std::vector<int>>& infected_node
        ,std::list<std::vector<int>>&  positions,std::mt19937& gen,std::uniform_int_distribution<int>& int_dist){
    int count(0),node;
    for (std::list< std::vector<int>>::iterator it_sigma=sigma.begin(),it_infected_node=infected_node.begin(),
                 it_positions=positions.begin();it_sigma!=sigma.end();++it_sigma,++it_infected_node,++it_positions){
        count=0;
        while (inital_infecteons>count){
//            node = rand()%N;
            node = int_dist(gen);
            if ((*it_sigma)[node]==0)
            {
                count++;
                (*it_sigma)[node]=1;
                (*it_positions)[node] = it_infected_node->size();
                it_infected_node->push_back(node);
            }
        }
    }
}

void intalize_infected_neighbours(int N,std::list< std::vector<int>> &infected,std::list< std::vector<int>> &infected_neighbors_in,
                                  std::list< std::vector<int>> &infected_neighbors_out,std::vector<int> &degrees_in
        ,std::vector<int>& degrees_out,std::vector<std::vector<int>>& Adjlist_in,std::vector<std::vector<int>>& Adjlist_out){
    for (std::list< std::vector<int>>::iterator it_infected_neighbors_in=infected_neighbors_in.begin(),
                 it_infected_neighbors_out=infected_neighbors_out.begin(),it_infected=infected.begin();it_infected_neighbors_in!=infected_neighbors_in.end()
            ;++it_infected_neighbors_in,++it_infected_neighbors_out,++it_infected){
        for (int i=0;i<it_infected->size();i++){
            increment_neighbours(degrees_out[(*it_infected)[i]],Adjlist_out[(*it_infected)[i]],*it_infected_neighbors_in);
            increment_neighbours(degrees_in[(*it_infected)[i]],Adjlist_in[(*it_infected)[i]],*it_infected_neighbors_out);
        }
    }
}


void intalized_sus(int N,std::list< std::vector<std::vector<int>>> &susceptible_nodes,std::list<std::vector<int>> &positions,std::list< std::vector<int>> &sigma,
                   std::list<std::vector<int>> &s_m,std::list< std::vector<int>> &infected_neighbors_in){
    std::list< std::vector<std::vector<int>>>::iterator it_susceptible_nodes=susceptible_nodes.begin();
    std::list<std::vector<int>>::iterator it_positions=positions.begin();
    std::list< std::vector<int>>::iterator it_sigma=sigma.begin();
    std::list<std::vector<int>>::iterator it_s_m=s_m.begin();
    std::list< std::vector<int>>::iterator it_infected_neighbors_in=infected_neighbors_in.begin();
    while (it_susceptible_nodes!=susceptible_nodes.end()){
        for (int i=0;i<N;i++){
            if ((*it_sigma)[i]==0){
                (*it_positions)[i] = (*it_susceptible_nodes)[(*it_infected_neighbors_in)[i]].size();
                (*it_susceptible_nodes)[(*it_infected_neighbors_in)[i]].push_back(i);
                (*it_s_m)[(*it_infected_neighbors_in)[i]]++;
            }
        }
        ++it_susceptible_nodes;
        ++it_positions;
        ++it_sigma;
        ++it_s_m;
        ++it_infected_neighbors_in;
    }
}

template <typename T>
void write_output_data(const std::string &filename,std::vector<T> &data){
    std::ofstream file(filename,std::ios::app);
    for (auto element=data.begin();element!=data.end();++element){file<<*element<<",";}
    file<<std::endl;
    file.close();
}

template <typename T>
void write_output_data(const std::string &filename,std::deque<T> &data){
    std::ofstream file(filename,std::ios::app);
    for (auto element=data.begin();element!=data.end();++element){file<<*element<<",";}
    file<<std::endl;
    file.close();
}

template <typename T1, typename T2>
void write_output_data_multi_col(const std::string &filename, const std::list<T1> &data1, std::list<T2> &data2){

    std::ofstream file(filename, std::ios::app);
    auto element1=data1.begin();
    auto element2=data2.begin();
    for (;element1!=data1.end() && element2!=data2.end();++element1,++element2){
        file<<*element1<<","<<*element2<<",";
        file<<std::endl;
    }
    file<<std::endl;
    file.close();
}

void intalize_SI_connections(std::list<int>& SI,std::list<std::vector<int>>& s_m,int k_max){
    std::list<std::vector<int>>::iterator it_s_m=s_m.begin();
    std::list<int>::iterator it_SI=SI.begin();
    while(it_SI!=SI.end()){
        for (int i=0;i<=k_max;i++) {
            (*it_SI) += (*it_s_m)[i] * i;
        }
        ++it_SI;
        ++it_s_m;
    }
}

void inital_networks_stat(int N,double x,int sims,std::list<int>& num_inf,std::list<std::vector<int>>& infected_node,
                          std::list< std::vector<int>>& infected_neighbors_in,std::list< std::vector<int>>& infected_neighbors_out,
                          std::list< std::vector<int>>& sigma,std::vector<int>& degrees_in,std::vector<int>& degrees_out,
                          std::vector<std::vector<int>>& Adjlist_in,std::vector<std::vector<int>>& Adjlist_out,std::list<std::vector<int>>&  positions
        ,std::list< std::vector<std::vector<int>>>& susceptible_nodes,std::list<int>& SI,std::list<std::vector<int>>& s_m,int k_max,
                          std::mt19937& gen,std::uniform_int_distribution<int>& int_dist){
    int inital_infecteons(int(x*N));
    num_inf = std::list<int>(sims,inital_infecteons);
    susceptible_nodes = std::list< std::vector<std::vector<int>>>(sims,std::vector<std::vector<int>>(k_max+1,std::vector<int>(0)));
    intalized_network_random_infection(N,inital_infecteons,sigma,infected_node,positions,gen,int_dist);
    intalize_infected_neighbours(N,infected_node,infected_neighbors_in,infected_neighbors_out,degrees_in,degrees_out,Adjlist_in,Adjlist_out);
    intalized_sus(N,susceptible_nodes,positions,sigma,s_m,infected_neighbors_in);
    intalize_SI_connections(SI,s_m,k_max);
}


//std::vector<double> GillespieMC(double tau,std::mt19937& gen,std::uniform_real_distribution<double>& uniform_dist,
//                   std::exponential_distribution<double> &exponential_dist,networks_dynamics &net_d,
//                   network_topology &net_t,double cat_start,double cat_duration){
//    simulation_data network(net_d);
//    std::pair<double,simulation_data> data_sim_gillespie(0.0,network);
//    std::vector<double> death;
//    while (!network.end(net_d)){
//        data_sim_gillespie = network.gillespie(tau,gen,uniform_dist,net_d,net_t,cat_start,cat_duration);
//        death.push_back(data_sim_gillespie.first);
//        if (data_sim_gillespie.first>0.0){continue;} // The simulation was erased so network points to the next network
//        ++network;
//    }
//    return death;
//}


std::pair<std::vector<double>,std::vector<double>> GillespieMC(double tau,std::mt19937& gen,std::uniform_real_distribution<double>& uniform_dist,
                   std::exponential_distribution<double> &exponential_dist,networks_dynamics &net_d,
                   network_topology &net_t,double cat_start,double cat_duration){
    simulation_data network(net_d);
    std::tuple<double, int, simulation_data> data_sim_gillespie(0.0,0,network);
    std::pair<double,int> sim_output;
    std::vector<double> extinction_time_num_inf,persistence_time_num_inf;
    int infected = 0;
    float time =0;
    while (!network.end(net_d)){
        data_sim_gillespie = network.gillespie(tau,gen,uniform_dist,net_d,net_t,cat_start,cat_duration);
        infected = std::get<1>(data_sim_gillespie);
        time = std::get<0>(data_sim_gillespie);
        if (infected>0){
            persistence_time_num_inf.push_back(time);
        }
        else{
            extinction_time_num_inf.push_back(time);
        }
        if (infected<=0){continue;} // The simulation was erased so network points to the next network
        ++network;
    }
    return std::make_pair(persistence_time_num_inf,extinction_time_num_inf);
}


std::pair<std::vector<double>,std::vector<double>> GillespieMC1d(double tau,std::mt19937& gen,std::uniform_real_distribution<double>& uniform_dist,
                                std::exponential_distribution<double> &exponential_dist,int N,
                                int seed,double cat_start,double cat_duration,int sim_max,double Alpha,double beta, double beta_cat,int k){
    std::vector<double> persistence,extinction;
    std::pair<double,int> network_extinction;
    for (int sim=0; sim<sim_max;sim++) {
        network_extinction = gillespie1d(tau, gen, uniform_dist, cat_start, cat_duration, N, seed, Alpha, beta,
                                         beta_cat, k);
        if (network_extinction.second > 0)
            persistence.push_back(network_extinction.first);
        else
            extinction.push_back(network_extinction.first);
//        death.push_back(gillespie1d(tau,gen,uniform_dist,cat_start,cat_duration,N, seed,Alpha, beta, beta_cat,k));
    }
    std::pair<std::vector<double>,std::vector<double>> output(persistence,extinction);
    return output;
}



std::deque<simulation_data>::iterator randomly_select_simulation(std::deque<simulation_data> &consmall
        ,std::deque<simulation_data>::iterator &end_simulation,double wtot,std::mt19937 &gen){
    std::uniform_real_distribution<> sim_dist(0.0,wtot);
    double weight_rand=sim_dist(gen),weight_sum=0.0;
    for (auto it_sim=consmall.begin();it_sim!=end_simulation;++it_sim){
        weight_sum+=*it_sim->i_weights;
        if (weight_sum>=weight_rand)
            return it_sim;
    }
    return end_simulation;
}


int main(int argc, char* argv[]) {
    std::string filename_in("/home/elad/multi_contact_rate_project/cpp_time_varying_rates/Adjin_0.txt"),
            filename_out("/home/elad/multi_contact_rate_project/cpp_time_varying_rates/Adjout_0.txt"),
            parametersname("/home/elad/multi_contact_rate_project/cpp_time_varying_rates/cparameters_0.txt");
    if (argc>1){
        filename_in=argv[1];
        filename_out=argv[2];
        parametersname = argv[3];
    }
    double death(0.0);
    std::fstream Adjfile_in,Adjfile_out,parametersfile;
    auto parameter_list=read_parameters(parametersname,parametersfile);
    int N=std::get<0>(parameter_list),sims=std::get<1>(parameter_list),start=std::get<2>(parameter_list),
            k=std::get<3>(parameter_list),network_number=std::get<8>(parameter_list);
    double inital_inf_percent=std::get<4>(parameter_list),lam=std::get<5>(parameter_list),Alpha=std::get<6>(parameter_list)
    ,Beta=std::get<7>(parameter_list),tau=std::get<9>(parameter_list),mf_solution=std::get<10>(parameter_list),
            cat_duration=std::get<16>(parameter_list),beta_cat=std::get<17>(parameter_list);
    std::string dir_path(std::get<13>(parameter_list)),prog=std::get<12>(parameter_list);
    int steps_c(sims),k_max(0),k_max_out(0); // Need to find ouf which dimension i np.size(n,1) refers to
    std::list<double> weights(sims,1.0/double(sims)),wg(sims,1.0/double(sims)),avec_sum(sims,0),t(sims,0);
    std::list<double>::iterator it_weights(weights.begin()),it_wg(wg.begin()),it_avec_sum(avec_sum.begin());
    std::vector<int> degrees_in(N,0),degrees_out(N,0);
    std::vector<std::vector<int>> Adjlist_in(N,std::vector<int>()),Adjlist_out(N,std::vector<int>());
    std::list<std::vector<int>> sigma(sims,std::vector<int>(N,0)),infected_node(sims,std::vector<int>()),
            ng(steps_c,std::vector<int>(N,0)),positions(sims,std::vector<int>(N,0));
    std::list<std::vector<int>>::iterator it_sigma(sigma.begin()),it_infected_node(infected_node.begin()),it_ng(ng.begin()),it_positions(positions.begin());
    std::list< std::vector<std::vector<int>>> susceptible_nodes(steps_c,std::vector<std::vector<int>>());
    std::list< std::vector<std::vector<int>>>::iterator it_susceptible_nodes(susceptible_nodes.begin());
    std::list<int> num_inf(sims,0),SI(sims,0);
    std::list<int>::iterator it_num_inf(num_inf.begin()),it_SI(SI.begin());
    std::list<double>::iterator it_t(t.begin());
    std::list< std::vector<int>> infected_neighbors_in(steps_c,std::vector<int>(N,0)),infected_neighbors_out(steps_c,std::vector<int>(N,0));
    std::list< std::vector<int>>::iterator it_infected_neighbors_in(infected_neighbors_in.begin()),it_infected_neighbors_out(infected_neighbors_out.begin());

    // Define distributions for generating random numbers
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    std::exponential_distribution<double> exponential_dist(1.0);
    std::uniform_int_distribution<int> int_dist(0, N - 1);


    // Ensure that there will be a new generation of random numbers in each iteration
    std::random_device rd;
    std::mt19937 gen(rd());

    auto start_time = std::chrono::high_resolution_clock::now();
    //Read the python's network structure, what are the different connections

    read_in_neighborslist(N,filename_in,Adjlist_in,degrees_in, Adjfile_in,k_max);
    read_in_neighborslist(N,filename_out,Adjlist_out,degrees_out,Adjfile_out,k_max_out);

    // Intailze the different networks list such as the number of infected nodes and so on
    std::list<std::vector<int>> s_m(steps_c,std::vector<int>(k_max+1,0));
    std::list<std::vector<int>>::iterator it_s_m(s_m.begin());
    inital_networks_stat(N,inital_inf_percent,sims,num_inf,infected_node,infected_neighbors_in,infected_neighbors_out,sigma,
                         degrees_in,degrees_out,Adjlist_in,Adjlist_out,positions,susceptible_nodes,SI,s_m,k_max,gen,int_dist);
//    std::vector<double> death_vec;
    std::pair<std::vector<double>,std::vector<double>> death_vec;
//    std::deque<double> Nlimits={0,mf_solution,N+1.0};

    std::list<int>::iterator it_min_new;
    int relaxation_time=20;
    // End of variable defnitions

    networks_dynamics net_d(num_inf,weights,avec_sum,t,infected_node,infected_neighbors_in,
                            infected_neighbors_out,sigma,s_m,positions,susceptible_nodes,SI);
    network_topology net_t(Alpha,Beta,k_max,degrees_in,degrees_out,Adjlist_in,Adjlist_out,beta_cat);

//    death_vec = GillespieMC(tau,gen,uniform_dist,exponential_dist,net_d,net_t,cat_start,cat_duration);
    if (eps_din!=0.0) {
        death_vec = GillespieMC(tau, gen, uniform_dist, exponential_dist, net_d, net_t, start, cat_duration);
    } else {
        death_vec = GillespieMC1d(tau, gen, uniform_dist, exponential_dist, N, int(inital_inf_percent * N), start,
                                  cat_duration, sims, Alpha, Beta, beta_cat, k);
    }
//    std::pair<std::vector<double>,std::vector<double>>  persistence_extinction = GillespieMC1d(tau,gen,uniform_dist,exponential_dist,N,int(inital_inf_percent*N),start,cat_duration,sims,Alpha,Beta,beta_cat,k);


//    # pragma omp parallel for //code for parallelization of jobs on cluster
//    for(int j=0;j<1;j++){
        // Run Gillespie's time step and update the number of deaths and net_d
        // death now contains the time of the simulations that died.
//        death_vec = GillespieMC(tau,gen,uniform_dist,exponential_dist,net_d,net_t,cat_start,cat_duration);
//        death_vec.push_back(death);

//        std::cout <<j<<std::endl;
//    }

//    auto start_pos_relax=death_vec.begin();
//    std::advance(start_pos_relax,relaxation_time);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    std::cout << "Program execution time: " << duration.count() << " seconds" << std::endl;
    std::string extinction_name="death.csv", parmetername="output_parameters.csv",persistence_name="persistence.csv";
    write_output_data(persistence_name,death_vec.first);
    write_output_data(extinction_name,death_vec.second);
//    write_output_data(name_Nlimits,Nlimits);
//    write_output_data_multi_col(weights_name,net_d.num_inf,net_d.weights);
    std::vector<double> outputparameters={double(duration.count()),double(network_number)};
    write_output_data(parmetername,outputparameters);
    return 0;
}