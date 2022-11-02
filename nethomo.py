import os
import networkx as nx
import numpy as np
import netinithomo
import rand_networks
import csv
import pickle

if __name__ == '__main__':
    Epsilon_sus = [0.0]
    Epsilon_inf = [0.0]
    # Epsilon=[0.02]
    eps_din,eps_dout = 0.0,0.0
    # eps_sus,eps_lam = 0.3,-0.3
    N = 2000
    k = 50
    x = 0.2
    Num_inf = int(x * N)
    Alpha = 1.0
    susceptibility = 'bimodal'
    infectability = 'bimodal'
    directed_model='uniform_c'
    prog = 'wnr' #can be either 'i' for the inatilization and reaching eq state or 'r' for running and recording fluc
    Lam = 1.5
    Time_limit = 400
    Start_recording_time = 100
    Beta_avg = Lam / k
    Num_different_networks= 10
    Num_inital_conditions= 100
    bank = 1000000
    parts = 1
    foldername ='track_inf_homo_nnn_N2000_k50_lam15_start100_end400_strength02'
    graphname  = 'GNull'
    count = 0
    susceptibility_avg = 1.0
    infectability_avg = 1.0
    sus_inf_correlation = 'ac'
    factor, duration, time_q,beta_time_type = 0.75, 10.0, 100.0,'c'
    rate_type ='c'

    if prog == 'i' or prog=='bi' or prog == 'si' or prog=='e' or prog=='ec' or prog=='ac' or prog=='r' or prog=='ri' or\
            prog=='g' or prog=='rg' or prog=='bd' or prog=='co' or prog=='cr' or prog=='q':
        os.mkdir(foldername)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(foldername)
    if prog=='i' or prog=='bi':
        G = nx.random_regular_graph(k,N)
        nx.write_gpickle(G,graphname)
    if prog == 'i':
        for epsilon in Epsilon:
            outfile ='i'+ str(epsilon).replace('.','')
            os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py ' + str(prog)+ ' ' +
                      str(epsilon) + ' ' + str(Beta_avg) +' ' + str(x) + ' ' + str(N) +
                      ' ' +str(graphname)+ ' ' + str(Alpha) + ' ' + str(Start_recording_time)+' '
                      + str(bank) + ' ' + str(outfile)+' '+str(type))
    elif prog=='r':
        for epsilon_sus,epsilon_inf in zip(Epsilon_sus,Epsilon_inf):
            Beta=Beta_avg/(1-epsilon_sus*epsilon_inf)
            for n in range(Num_different_networks):
                beta_inf,beta_sus=netinithomo.bi_beta_anti_correlated(N,epsilon_inf,epsilon_sus,1.0)
                G = nx.random_regular_graph(k, N)
                # G = nx.complete_graph(N)
                G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
                infile = graphname + '_' + str(epsilon_sus).replace('.', '') + '_' + str(n)+'.pickle'
                nx.write_gpickle(G, infile)
                outfile ='o'+str(epsilon_sus).replace('.', '')
                for p in range(parts):
                    os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                              str(Alpha) + ' ' + str(Time_limit)+ ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))

    elif prog=='ri':
        for epsilon_sus,epsilon_inf in zip(Epsilon_sus,Epsilon_inf):
            Beta=Beta_avg/(1-epsilon_sus*epsilon_inf)
            for n in range(Num_different_networks):
                beta_inf,beta_sus=netinithomo.bi_beta_anti_correlated(N,epsilon_inf,epsilon_sus,1.0)
                # G = nx.random_regular_graph(k, N)
                G = nx.complete_graph(N)
                G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
                infile = graphname + '_' + str(epsilon_sus).replace('.', '') + '_' + str(n)+'.pickle'
                nx.write_gpickle(G, infile)
                outfile ='o'+str(epsilon_sus).replace('.', '')
                for p in range(parts):
                    os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                              str(Alpha) + ' ' + str(Time_limit)+ ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta) + ' '+str(Start_recording_time))
    elif prog=='bi':
        for epsilon in Epsilon:
            outfile = 'i' + str(epsilon).replace('.', '')
            os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py ' + str(prog) + ' ' +
                      str(epsilon) + ' ' + str(Beta_avg) + ' ' + str(x) + ' ' + str(N) +
                      ' ' + str(graphname) + ' ' + str(Alpha) + ' ' + str(Start_recording_time) + ' '
                      + str(bank) + ' ' + str(outfile))
    elif prog=='b':
        for epsilon in Epsilon:
            infile = 'i'+str(epsilon).replace('.', '')
            outfile ='o'+str(epsilon).replace('.', '')
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(Time_limit) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile)+
                          ' '+str(epsilon)+ ' '+ str(Beta_avg))
    elif prog=='si':
        for epsilon in Epsilon:
            count=count+1
            G =nx.random_regular_graph(k,N)
            nx.write_gpickle(G,graphname+'_'+str(count))
            outfile ='i'+ str(epsilon).replace('.','')
            os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py ' + str(prog)+ ' ' +
                      str(epsilon) + ' ' + str(Beta_avg) +' ' + str(x) + ' ' + str(N) +
                      ' ' +str(graphname)+ ' ' + str(Alpha) + ' ' + str(Start_recording_time)+' '
                      + str(bank) + ' ' + str(outfile)+' '+str(type)+ ' ' + '_'+str(count))
    elif prog =='sr':
        for epsilon in Epsilon:
            count  = count +1
            infile = 'i'+str(epsilon).replace('.', '')
            outfile ='o'+str(epsilon).replace('.', '')
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(Time_limit) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' +'_'+str(count))
    elif prog =='e':
        for epsilon_sus,epsilon_inf in zip(Epsilon_sus,Epsilon_inf):
            Beta=Beta_avg/(1+epsilon_sus*epsilon_inf)
            for n in range(Num_different_networks):
                G = nx.random_regular_graph(k, N)
                choose_beta = lambda net_dist,avg,epsilon: np.random.normal(avg, epsilon * avg, N) \
                    if net_dist == 'gauss' else np.random.gamma((avg / epsilon) ** 2, epsilon ** 2 / avg, N) \
                    if net_dist == 'gamma' else np.zeros(N) if net_dist =='z' else np.ones(N) if net_dist =='ones' else netinithomo.bi_beta(N, epsilon, avg)
                beta_sus = choose_beta(susceptibility,susceptibility_avg,epsilon_sus)
                beta_inf = choose_beta(infectability,infectability_avg,epsilon_inf)
                # beta_all = choose_beta(susceptibility,susceptibility_avg,epsilon_sus)
                G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
                infile = graphname + '_' + str(epsilon_sus).replace('.', '') + '_' + str(n)+'.pickle'
                nx.write_gpickle(G, infile)
                outfile ='o'+str(epsilon_sus).replace('.', '')
                for p in range(parts):
                    os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                              str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='ec':
        for epsilon_sus,epsilon_inf in zip(Epsilon_sus,Epsilon_inf):
            Beta=Beta_avg/(1+epsilon_sus*epsilon_inf)
            for n in range(Num_different_networks):
                beta_inf,beta_sus=netinithomo.bi_beta_correlated(N,epsilon_inf,epsilon_sus,1.0)
                # G = nx.random_regular_graph(k, N)
                G = nx.complete_graph(N)
                G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
                infile = graphname + '_' + str(epsilon_sus).replace('.', '') + '_' + str(n)+'.pickle'
                nx.write_gpickle(G, infile)
                outfile ='o'+str(epsilon_sus).replace('.', '')
                for p in range(parts):
                    os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                              str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='ac':
        for epsilon_sus,epsilon_inf in zip(Epsilon_sus,Epsilon_inf):
            Beta=Beta_avg/(1-epsilon_sus*epsilon_inf)
            for n in range(Num_different_networks):
                beta_inf,beta_sus=netinithomo.bi_beta_anti_correlated(N,epsilon_inf,epsilon_sus,1.0)
                # G = nx.random_regular_graph(k, N)
                G = nx.complete_graph(N)
                G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
                infile = graphname + '_' + str(epsilon_sus).replace('.', '') + '_' + str(n)+'.pickle'
                nx.write_gpickle(G, infile)
                outfile ='o'+str(epsilon_sus).replace('.', '')
                for p in range(parts):
                    os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                              str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='g':
        # Program to find extinction time in case I have genreal rate distubtion (tringular etc)
        for epsilon_sus,epsilon_inf in zip(Epsilon_sus,Epsilon_inf):
            Beta=Beta_avg/(1-epsilon_sus*epsilon_inf) if sus_inf_correlation == 'a' else Beta_avg/(1+epsilon_sus*epsilon_inf)
            for n in range(Num_different_networks):
                # beta_inf,beta_sus=netinithomo.triangular_beta(N,epsilon_inf,epsilon_sus,sus_inf_correlation)
                beta_inf,beta_sus=netinithomo.uniform_beta(N,epsilon_inf,epsilon_sus,sus_inf_correlation)
                # G = nx.random_regular_graph(k, N)
                G = nx.complete_graph(N)
                G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
                infile = graphname + '_' + str(epsilon_sus).replace('.', '') + '_' + str(n)+'.pickle'
                nx.write_gpickle(G, infile)
                outfile ='o'+str(epsilon_sus).replace('.', '')
                for p in range(parts):
                    os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                              str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog=='rg':
        for epsilon_sus,epsilon_inf in zip(Epsilon_sus,Epsilon_inf):
            Beta=Beta_avg/(1-epsilon_sus*epsilon_inf) if sus_inf_correlation == 'a' else Beta_avg/(1+epsilon_sus*epsilon_inf)
            for n in range(Num_different_networks):
                # beta_inf,beta_sus=netinithomo.triangular_beta(N,epsilon_inf,epsilon_sus,sus_inf_correlation)
                beta_inf,beta_sus=netinithomo.uniform_beta(N,epsilon_inf,epsilon_sus,sus_inf_correlation)
                # G = nx.random_regular_graph(k, N)
                G = nx.complete_graph(N)
                G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
                infile = graphname + '_' + str(epsilon_sus).replace('.', '') + '_' + str(n)+'.pickle'
                nx.write_gpickle(G, infile)
                outfile ='o'+str(epsilon_sus).replace('.', '')
                for p in range(parts):
                    os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                              str(Alpha) + ' ' + str(Time_limit)+ ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta) + ' '+str(Start_recording_time))
    elif prog =='bd':
        Beta=Beta_avg/(1+eps_din*eps_dout)
        d1_in, d1_out, d2_in, d2_out =int(k*(1-eps_din)),int(k*(1-eps_dout)),int(k*(1+eps_din)),int(k*(1+eps_dout))
        for n in range(Num_different_networks):
            G = rand_networks.random_bimodal_directed_graph(d1_in, d1_out,d2_in,d2_out,N)
            G = netinithomo.set_graph_attriubute_DiGraph(G)
            infile = graphname + '_' + str(eps_din).replace('.', '') + '_' + str(n)+'.pickle'
            nx.write_gpickle(G, infile)
            outfile ='o_d1in' + str(d1_in).replace('.', '') +'_o_d1out' + str(d1_out).replace('.', '')
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='co':
        for n in range(Num_different_networks):
            G = rand_networks.configuration_model_directed_graph(directed_model, eps_din,eps_dout,k,N)
            G = netinithomo.set_graph_attriubute_DiGraph(G)
            infile = graphname + '_' + str(eps_din).replace('.', '') + '_' + str(n)+'.pickle'
            nx.write_gpickle(G, infile)
            outfile ='o_eps_in' + str(np.abs(eps_din)).replace('.', '') +'eps_dout' + str(np.abs(eps_dout)).replace('.', '')
            k_avg_graph = np.mean([G.in_degree(n) for n in G.nodes()])
            Beta_graph = Lam/k_avg_graph
            eps_in_graph = np.std([G.in_degree(n) for n in G.nodes()])/k_avg_graph
            eps_out_graph = np.std([G.out_degree(n) for n in G.nodes()])/k_avg_graph
            Beta = Beta_graph / (1 + np.sign(eps_din)*eps_in_graph * np.sign(eps_dout)* eps_out_graph)
            f = open('parameters_'+outfile + '.csv', "a+")
            with f:
                writer = csv.writer(f)
                writer.writerows([[k_avg_graph, np.sign(eps_din)*eps_in_graph,np.sign(eps_dout)*eps_out_graph]])
            f.close()
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='cr':
        for n in range(Num_different_networks):
            beta_inf,beta_sus=netinithomo.general_beta(N,eps_lam,eps_sus,directed_model,k)
            G = nx.random_regular_graph(k, N)
            # G = nx.complete_graph(N)
            G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
            infile = graphname + '_' + str(eps_sus).replace('.', '') + '_' + str(n)+'.pickle'
            nx.write_gpickle(G, infile)
            eps_sus_graph = np.std(beta_sus)/np.mean(beta_sus)
            eps_lam_graph = np.std(beta_inf)/np.mean(beta_inf)
            Beta = Beta_avg / (1 + np.sign(eps_sus)*eps_sus_graph * np.sign(eps_lam)* eps_lam_graph)
            outfile ='o'+str(eps_sus).replace('.', '')
            f = open('parameters_'+outfile + '.csv', "a+")
            with f:
                writer = csv.writer(f)
                writer.writerows([[np.mean(beta_sus), np.sign(eps_sus)*eps_sus_graph,np.sign(eps_lam)*eps_lam_graph]])
            f.close()
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='q':
        # run until either time t is reached or extinction, this program record how many extinction events occured during that time in adaptive network
        Beta = Beta_avg / (1 + eps_din * eps_dout)
        d1_in, d1_out, d2_in, d2_out = int(k * (1 - eps_din)), int(k * (1 - eps_dout)), int(k * (1 + eps_din)), int(
            k * (1 + eps_dout))
        for n in range(Num_different_networks):
            G = rand_networks.random_bimodal_directed_graph(d1_in, d1_out,d2_in,d2_out,N)
            G = netinithomo.set_graph_attriubute_DiGraph(G)
            # beta_inf,beta_sus = netinithomo.bi_beta_correlated(N,0.0,0.0,1.0)
            # G = nx.random_regular_graph(k, N)
            # G = nx.complete_graph(N)
            infile = graphname + '_' + str(eps_din).replace('.', '') + '_' + str(n)+'.pickle'
            nx.write_gpickle(G, infile)
            outfile ='o_d1in' + str(d1_in).replace('.', '') +'_o_d1out' + str(d1_out).replace('.', '')
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(Time_limit) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' +
                          str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta) +
                          ' ' + str(factor) + ' ' + str(duration) + ' ' + str(time_q) + ' ' + str(beta_time_type))
    elif prog =='th':
        beta=lambda t: Beta_avg
        for n in range(Num_different_networks):
            G = nx.random_regular_graph(k, N)
            # G = nx.complete_graph(N)
            G = netinithomo.intalize_homo_temporal_graph(G, N, beta)
            infile = graphname + '_' + str(Beta_avg).replace('.', '') + '_' + str(n)+'.pickle'
            nx.write_gpickle(G, infile)
            outfile ='o'+str(Beta_avg).replace('.', '')
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(Time_limit) + ' ' + str(bank)+ ' ' + str(outfile) + ' ' +
                          str(infile) + ' ' + str(Num_inital_conditions)
                          + ' ' +str(Num_inf) + ' ' +str(n) + ' ' +str(Start_recording_time)+' ' +str(rate_type))