import numpy as np
import bisect
import netinithomo
import networkx as nx
import csv
import pickle
import sys
import rand_networks

import matplotlib.pyplot as plt



def fluctuation_init(epsilon,avg_beta,x,N,G_name,Alpha,Time_limit,bank,outfilename,type,ending=''):
    choose_beta = lambda type: np.random.normal(avg_beta, epsilon * avg_beta, N) \
        if type == 'gauss' else np.random.gamma((avg_beta / epsilon) ** 2, epsilon ** 2 / avg_beta, N) \
        if type == 'gamma' else netinithomo.bi_beta(N, epsilon, avg_beta)
    Beta = choose_beta(type)
    G = nx.read_gpickle(G_name+ending)
    T = []
    Num_inf = int(x * N)
    R_tot, Rates = netinithomo.intialize_graph(G, N, Num_inf, Beta, Alpha)
    Total_time = 0.0
    T.append(Total_time)

    count = 0
    r = np.random.uniform(0, 1, (bank, 2))
    ######################
    # Main Gillespie Loop
    ######################
    while Num_inf > 0 and Total_time<Time_limit:
        R_norm = np.cumsum(Rates)
        r_pos = R_tot * r[count, 1]
        person = bisect.bisect_left(R_norm, r_pos)

        tau= np.log(1 / r[count, 0]) / R_tot
        Total_time = Total_time + tau

        # rand_networks.draw_infections_nx_g(G,pos_nodes_plt,'frame'+str(count),[person])

        if G.nodes[person]['infected'] == True:
            Num_inf = Num_inf - 1
            Rates[person] = 0.0
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] - G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot - G.nodes[Neighbor]['contact_rate']
                else:
                    Rates[person] = Rates[person] + G.nodes[person]['contact_rate']
            R_tot = R_tot + Rates[person] - Alpha
            G.nodes[person]['infected'] = False
        else:
            Num_inf = Num_inf + 1
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] + G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot + G.nodes[Neighbor]['contact_rate']
            R_tot = R_tot - Rates[person] + Alpha
            Rates[person] = Alpha
            G.nodes[person]['infected'] = True
        count = count + 1
        if count >= bank:
            r = np.random.uniform(0, 1, (bank, 2))
            count = 0
    outfile = open (outfilename + ending + '_Rates.pickle','wb')
    pickle.dump(Rates,outfile)
    outfile.close()
    outfile = open(outfilename + ending + '_R_tot.pickle','wb')
    pickle.dump(R_tot, outfile)
    outfile.close()
    outfile = open(outfilename + ending + '_Num_inf.pickle','wb')
    pickle.dump(Num_inf, outfile)
    outfile.close()
    nx.write_gpickle(G,outfilename + ending + '_G.pickle')
    return 0

def fluctuation_run(Alpha,Time_limit,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G=nx.read_gpickle(infile)
    seed_nodes= Num_inf
    for run_loop_counter in range(runs):
        T = []
        I = []
        runs_csv=[]
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)

        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)

        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                    else:
                        Rates[person] = Rates[person] + Beta*(G.nodes[person]['contact_rate'] * G.nodes[Neighbor]['spread_rate'])
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
            if Total_time-T[-1]>=0.1:
                I.append(Num_inf)
                T.append(Total_time)
                net_num.append(network_number)
                runs_csv.append(run_loop_counter)
        f = open(outfile + '.csv', "a+")
        l = [T, I, runs_csv,net_num]
        l = zip(*l)
        with f:
            writer = csv.writer(f)
            writer.writerows(l)
        f.close()
    return 0


def fluctuation_run_catastrophe(Alpha,Time_limit,bank,outfile,infile,runs,Num_inf,network_number,Beta,factor,duration,time_q,type_q):

    def create_time_varying_beta(Total_time,time_q,duration,quarntine,G,Alpha,Beta,type_q,R_tot,Rates,Beta_org):
        if type_q=='c':
            if Total_time > time_q and Total_time <= time_q + duration and quarntine == False:
                # start a quarntine if there isn't one already and if the time is right
                Beta = Beta * factor
                R_tot, Rates = netinithomo.inatlize_quarntine_graph(G, G.number_of_nodes(), Alpha, Beta)
                quarntine = True
            elif Total_time > time_q and Total_time > time_q + duration and quarntine == True:
                # stop the quarantine and resume the previous infection rate
                Beta = Beta_org
                R_tot, Rates = netinithomo.inatlize_quarntine_graph(G, G.number_of_nodes(), Alpha, Beta)
                quarntine = False
        elif type_q == 'p':
            if Total_time > time_q and Total_time <= time_q + duration:
                Beta = Beta * (np.cos((2 * np.pi * (Total_time - time_q)) / factor)) ** 2
                R_tot, Rates = netinithomo.inatlize_quarntine_graph(G, G.number_of_nodes(), Alpha, Beta)
            elif Total_time > time_q and Total_time > time_q + duration:
                # stop the quarantine and resume the previous infection rate
                Beta = Beta_org
                R_tot, Rates = netinithomo.inatlize_quarntine_graph(G, G.number_of_nodes(), Alpha, Beta)
        return quarntine,Beta, R_tot, Rates

    G=nx.read_gpickle(infile)
    seed_nodes= Num_inf
    Beta_org=Beta
    for run_loop_counter in range(runs):
        T = []
        I = []
        runs_csv=[]
        Beta=Beta_org
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        nx.set_node_attributes(G,False,'infected')
        R_tot, Rates = netinithomo.inatlize_inf_DiGraph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)
        quarntine = False

        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit:

            quarntine,Beta, R_tot, Rates = create_time_varying_beta(Total_time,time_q,duration,quarntine,G,Alpha,Beta,type_q,R_tot,Rates,Beta_org)

            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta
                        R_tot = R_tot - Beta
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        Rates[person] = Rates[person] + Beta
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta
                        R_tot = R_tot + Beta
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number,Num_inf]])
        f.close()
    return 0


def fluctuation_run_no_decay(Alpha,Time_limit,bank,outfile,infile,runs,Num_inf,network_number,Beta,start_recording_time):
    G=nx.read_gpickle(infile)
    seed_nodes= Num_inf
    for run_loop_counter in range(runs):
        T = []
        I = []
        I_type_above_avg,I_type_below_avg= [],[]
        runs_csv=[]
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))

        for l in range(G.number_of_nodes()):
            G.nodes[l]['infected'] = False
        R_tot, Rates = netinithomo.inatlize_inf_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)

        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)

        num_inf_above_avg, num_inf_below_avg = 0, 0
        for node_type_number in range(G.number_of_nodes()):
            if G.nodes[node_type_number]['infected'] == True:
                if G.nodes[node_type_number]['contact_rate'] > 1.0:
                    num_inf_above_avg = num_inf_above_avg + 1
                else:
                    num_inf_below_avg = num_inf_below_avg + 1
        I_type_above_avg.append(num_inf_above_avg)
        I_type_below_avg.append(num_inf_below_avg)

        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True and Num_inf>1:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                    else:
                        Rates[person] = Rates[person] + Beta*(G.nodes[person]['contact_rate'] * G.nodes[Neighbor]['spread_rate'])
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            elif G.nodes[person]['infected'] == False :
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
            if Total_time-T[-1]>=0.1 and Total_time>=start_recording_time:
                I.append(Num_inf)
                T.append(Total_time)
                net_num.append(network_number)
                runs_csv.append(run_loop_counter)
                num_inf_above_avg,num_inf_below_avg=0,0
                for node_type_number in range(G.number_of_nodes()):
                    if G.nodes[node_type_number]['infected'] == True:
                        if G.nodes[node_type_number]['contact_rate']>1.0:
                            num_inf_above_avg = num_inf_above_avg + 1
                        else:
                            num_inf_below_avg = num_inf_below_avg + 1
                I_type_above_avg.append(num_inf_above_avg)
                I_type_below_avg.append(num_inf_below_avg)
        f = open(outfile + '.csv', "a+")
        l = [T, I,I_type_above_avg,I_type_below_avg ,runs_csv,net_num]
        l = zip(*l)
        with f:
            writer = csv.writer(f)
            writer.writerows(l)
        f.close()
    return 0


def fluctuation_run_extinction(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1


            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                    else:
                        Rates[person] = Rates[person] + Beta*(G.nodes[person]['contact_rate'] * G.nodes[Neighbor]['spread_rate'])
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
    return 0



def fluctuation_run_extinction_undirected_graph(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_undirected(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1


            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta
                        R_tot = R_tot - Beta
                    else:
                        Rates[person] = Rates[person] + Beta
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta
                        R_tot = R_tot + Beta
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
    return 0


def fluctuation_run_extinction_DiGraph(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_DiGraph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.node[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta
                        R_tot = R_tot - Beta
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        Rates[person] = Rates[person] + Beta
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta
                        R_tot = R_tot + Beta
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
    return 0


def fluctuation_init_track_lam(epsilon,avg_beta,x,N,G_name,Alpha,Time_limit,bank,outfilename):
    lam_plus, lam_minus = avg_beta*(1+epsilon), avg_beta*(1-epsilon)
    Beta = netinithomo.bi_beta(N, epsilon, avg_beta)
    G=nx.read_gpickle(G_name)
    T = []
    Num_inf = int(x * N)
    R_tot, Rates,Num_inf_plus, Num_inf_minus = netinithomo.intialize_graph_track_bimodal(G, N, Num_inf, Beta, Alpha,lam_plus,lam_minus)
    Total_time = 0.0
    T.append(Total_time)

    count = 0
    r = np.random.uniform(0, 1, (bank, 2))
    ######################
    # Main Gillespie Loop
    ######################
    while Num_inf > 0 and Total_time<Time_limit:
        R_norm = np.cumsum(Rates)
        r_pos = R_tot * r[count, 1]
        person = bisect.bisect_left(R_norm, r_pos)

        tau= np.log(1 / r[count, 0]) / R_tot
        Total_time = Total_time + tau

        # rand_networks.draw_infections_nx_g(G,pos_nodes_plt,'frame'+str(count),[person])

        if G.nodes[person]['infected'] == True:
            Num_inf = Num_inf - 1
            if round(G.nodes[person]['contact_rate'],6)==round(lam_plus,6):
                Num_inf_plus = Num_inf_plus -1
            else:
                Num_inf_minus = Num_inf_minus - 1
            Rates[person] = 0.0
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] - G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot - G.nodes[Neighbor]['contact_rate']
                else:
                    Rates[person] = Rates[person] + G.nodes[person]['contact_rate']
            R_tot = R_tot + Rates[person] - Alpha
            G.nodes[person]['infected'] = False
        else:
            Num_inf = Num_inf + 1
            if round(G.nodes[person]['contact_rate'],6)==round(lam_plus,6):
                Num_inf_plus = Num_inf_plus + 1
            else:
                Num_inf_minus = Num_inf_minus + 1
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] + G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot + G.nodes[Neighbor]['contact_rate']
            R_tot = R_tot - Rates[person] + Alpha
            Rates[person] = Alpha
            G.nodes[person]['infected'] = True
        count = count + 1
        if count >= bank:
            r = np.random.uniform(0, 1, (bank, 2))
            count = 0
    outfile = open (outfilename+'_Rates.pickle','wb')
    pickle.dump(Rates,outfile)
    outfile.close()
    outfile = open(outfilename + '_R_tot.pickle','wb')
    pickle.dump(R_tot, outfile)
    outfile.close()
    outfile = open(outfilename+'_Num_inf.pickle','wb')
    pickle.dump(Num_inf, outfile)
    outfile.close()
    outfile = open(outfilename+'_Num_inf_plus.pickle','wb')
    pickle.dump(Num_inf_plus, outfile)
    outfile.close()
    outfile = open(outfilename+'_Num_inf_minus.pickle','wb')
    pickle.dump(Num_inf_minus, outfile)
    outfile.close()
    nx.write_gpickle(G,outfilename+'_G.pickle')
    return 0


def fluctuation_run_track_lam(Alpha,Time_limit,bank,outfile,infile,epsilon,avg_beta):
    lam_plus, lam_minus = avg_beta*(1+epsilon), avg_beta*(1-epsilon)
    G=nx.read_gpickle(infile+'_G.pickle')
    with open(infile+'_Num_inf.pickle','rb') as pickle_file:
        Num_inf = pickle.load(pickle_file)
    pickle_file.close()
    with open(infile+'_Num_inf_plus.pickle','rb') as pickle_file:
        Num_inf_plus = pickle.load(pickle_file)
    pickle_file.close()
    with open(infile+'_Num_inf_minus.pickle','rb') as pickle_file:
        Num_inf_minus = pickle.load(pickle_file)
    pickle_file.close()
    with open(infile + '_R_tot.pickle','rb') as pickle_file:
        R_tot = pickle.load(pickle_file)
    pickle_file.close()
    with open(infile + '_Rates.pickle','rb') as pickle_file:
        Rates = pickle.load(pickle_file)
    pickle_file.close()
    T, I, I_plus, I_minus = [], [], [], []
    I.append(Num_inf)
    I_plus.append(Num_inf_plus)
    I_minus.append(Num_inf_minus)
    Total_time = 0.0
    T.append(Total_time)

    count = 0
    r = np.random.uniform(0, 1, (bank, 2))
    ######################
    # Main Gillespie Loop
    ######################
    while Num_inf > 0 and Total_time<Time_limit:
        R_norm = np.cumsum(Rates)
        r_pos = R_tot * r[count, 1]
        person = bisect.bisect_left(R_norm, r_pos)

        tau= np.log(1 / r[count, 0]) / R_tot
        Total_time = Total_time + tau

        # rand_networks.draw_infections_nx_g(G,pos_nodes_plt,'frame'+str(count),[person])

        if G.nodes[person]['infected'] == True:
            Num_inf = Num_inf - 1
            if round(G.nodes[person]['contact_rate'],6)==round(lam_plus,6):
                Num_inf_plus = Num_inf_plus - 1
            else:
                Num_inf_minus = Num_inf_minus -1
            Rates[person] = 0.0
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] - G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot - G.nodes[Neighbor]['contact_rate']
                else:
                    Rates[person] = Rates[person] + G.nodes[person]['contact_rate']
            R_tot = R_tot + Rates[person] - Alpha
            G.nodes[person]['infected'] = False
        else:
            Num_inf = Num_inf + 1
            if round(G.nodes[person]['contact_rate'],6)==round(lam_plus,6):
                Num_inf_plus = Num_inf_plus + 1
            else:
                Num_inf_minus = Num_inf_minus +1
            for Neighbor in G[person]:
                if G.nodes[Neighbor]['infected'] == False:
                    Rates[Neighbor] = Rates[Neighbor] + G.nodes[Neighbor]['contact_rate']
                    R_tot = R_tot + G.nodes[Neighbor]['contact_rate']
            R_tot = R_tot - Rates[person] + Alpha
            Rates[person] = Alpha
            G.nodes[person]['infected'] = True
        count = count + 1
        if count >= bank:
            r = np.random.uniform(0, 1, (bank, 2))
            count = 0
        if Total_time-T[-1]>=1:
            I.append(Num_inf)
            I_plus.append(Num_inf_plus)
            I_minus.append(Num_inf_minus)
            T.append(round(Total_time))
    f = open(outfile+'_total.csv',"a+")
    with f:
        writer = csv.writer(f)
        writer.writerows([I])
    f.close()
    f = open(outfile+'infected_plus.csv',"a+")
    with f:
        writer = csv.writer(f)
        writer.writerows([I_plus])
    f.close()
    f = open(outfile+'infected_minus.csv',"a+")
    with f:
        writer = csv.writer(f)
        writer.writerows([I_minus])
    f.close()
    return 0


def actasmain():
    Epsilon_sus = [0.0]
    Epsilon_inf = [0.0]
    Epsilon=[0.0]
    N = 100
    k = 50
    x = 0.2
    eps_din,eps_dout = 0.1,0.1
    eps_sus,eps_lam = 0.1,0.1
    Num_inf = int(x * N)
    Alpha = 1.0
    susceptibility = 'bimodal'
    infectability = 'bimodal'
    directed_model='gauss_c'
    prog = 'q' #can be either 'i' for the inatilization and reaching eq state or 'r' for running and recording fluc
    Lam = 1.6
    Time_limit = 200
    Start_recording_time = 50
    Beta_avg = Lam / k
    Num_different_networks= 2
    Num_inital_conditions= 2
    bank = 1000000
    parts = 1
    graphname  = 'GNull'
    count = 0
    susceptibility_avg = 1.0
    infectability_avg = 1.0
    foldername ='base'
    graphname  = 'GNull'
    outfile ='o'
    sus_inf_correlation = 'a'
    # Beta = Beta_avg / (1 + Epsilon_sus[0] * Epsilon_inf[0])
    # Beta = Beta_avg / (1 - Epsilon_sus[0] * Epsilon_inf[0]) if sus_inf_correlation is 'a' else Beta_avg / (
    #             1 + Epsilon_sus[0] * Epsilon_inf[0])
    Beta = Beta_avg / (1 + eps_lam * eps_sus)
    factor, duration, time_q = 0.5, 0.5, 100.0
    beta_time_type='p'

    # G = nx.random_regular_graph(k, N)
    # G = nx.complete_graph(N)
    # beta_inf, beta_sus = netinithomo.general_beta(N, eps_lam, eps_sus, directed_model, k)
    # beta_inf, beta_sus = netinithomo.bi_beta_correlated(N, 0.0, 0.0, 1.0)
    # G = netinithomo.intalize_lam_graph(G, N, beta_sus, beta_inf)
    d1_in, d1_out, d2_in, d2_out = int(k * (1 - eps_din)), int(k * (1 - eps_dout)), int(k * (1 + eps_din)), int(
        k * (1 + eps_dout))
    G = rand_networks.random_bimodal_directed_graph(d1_in, d1_out, d2_in, d2_out, N)
    G = netinithomo.set_graph_attriubute_DiGraph(G)

    # choose_beta = lambda net_dist, avg, epsilon: np.random.normal(avg, epsilon * avg, N) \
    #     if net_dist == 'gauss' else np.random.gamma((avg / epsilon) ** 2, epsilon ** 2 / avg, N) \
    #     if net_dist == 'gamma' else np.zeros(N) if net_dist == 'z' else np.ones(
    #     N) if net_dist == 'ones' else netinithomo.bi_beta(N, epsilon, avg)

    # beta_sus = choose_beta(susceptibility, susceptibility_avg,Epsilon_sus[0])
    # beta_inf = choose_beta(infectability, infectability_avg,Epsilon_inf[0])
    # G = netinithomo.intalize_lam_graph(G, N, beta_sus, beta_inf)
    # beta_inf, beta_sus = netinithomo.triangular_beta(N, Epsilon_inf[0], Epsilon_sus[0], sus_inf_correlation)
    # G = rand_networks.configuration_model_directed_graph(directed_model,eps_dout,eps_din,k,N)
    # G = netinithomo.set_graph_attriubute_DiGraph(G)
    # G = rand_networks.configuration_model_undirected_graph(directed_model,eps_dout,eps_din,k,N)
    # G = netinithomo.set_graph_attriubute_DiGraph(G)
    # G = rand_networks.configuration_model_undirected_graph(0.1,k,N)
    # G = rand_networks.jason_graph('jason_trans_file_no_degree.csv')
    # G = netinithomo.set_graph_attriubute_DiGraph(G)
    # G = netinithomo.intalize_lam_graph(G, N, np.ones(N), np.ones(N))
    # Beta = Beta_avg / (1 + eps_din * eps_dout)
    # Beta = Beta/np.mean([G.in_degree(n) for n in G.nodes()])
    # Beta = Lam/np.mean([G.in_degree(n) for n in G.nodes()])
    n=0


    # fluctuation_run_extinction_weighted_graph(Alpha, bank, outfile, infile,
    #                                           Num_inital_conditions, Num_inf, n,
    #                                           Beta)
    nx.write_gpickle(G, graphname)
    # infile = graphname + '_' + str(epsilon).replace('.', '') + '_' + str(n)+'.pickle'
    infile=graphname
    # fluctuation_run_extinction(Alpha,bank,outfile,infile,Num_inital_conditions,Num_inf,1,Beta)
    fluctuation_run_catastrophe(Alpha,Time_limit,bank,outfile,infile,Num_inital_conditions,Num_inf,n,Beta,factor,duration,time_q,beta_time_type)
    # fluctuation_run_no_decay(Alpha, Time_limit, bank, outfile, infile, Num_inital_conditions,
    #                 Num_inf, 1, Beta,Start_recording_time)
    # fluctuation_run_extinction_undirected_graph(Alpha, bank, outfile, infile, Num_inital_conditions,
    #                            Num_inf, 1, Beta)
    # fluctuation_run_extinction_DiGraph(Alpha, bank, outfile, infile, Num_inital_conditions,
    #                                    Num_inf, 1, Beta)

if __name__ == '__main__':
    submit = True
    if submit==False:
        actasmain()
    else:
         if sys.argv[1]=='i':
             # Parmeters order: epsilon,avg_beta,x,N,G_name,Alpha,Time_limit,bank,outfilename,type
             fluctuation_init(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), str(sys.argv[6]), float(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9]),sys.argv[10],sys.argv[11])
         elif sys.argv[1]=='r':
             # Parameters order: Alpha,Time_limit,bank,outfile,infile
             fluctuation_run(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6],int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]),float(sys.argv[10]))
         elif sys.argv[1]=='b':
             #Parameters: Alpha,Time_limit,bank,outfile,infile,epsilon,avg_beta
             fluctuation_run_track_lam(float(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],sys.argv[6],float(sys.argv[7]),float(sys.argv[8]))
         elif sys.argv[1]=='bi':
             fluctuation_init_track_lam(float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),int(sys.argv[5]),sys.argv[6],float(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]),sys.argv[10])
         elif sys.argv[1]=='si':
             # Parmeters order: epsilon,avg_beta,x,N,G_name,Alpha,Time_limit,bank,outfilename,type,ending
             fluctuation_init(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), str(sys.argv[6]), float(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9]),sys.argv[10],sys.argv[11],sys.argv[12])
         elif sys.argv[1]=='sr':
             # Parameters order: Alpha,Time_limit,bank,outfile,infile,ending
             fluctuation_run(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6],sys.argv[7])
         elif sys.argv[1]=='e' or sys.argv[1]=='ec' or sys.argv[1]=='ac' or sys.argv[1] == 'g' or sys.argv[1] == 'cr' :
             # Parameters order: Alpha,bank,outfile,infile
             fluctuation_run_extinction(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),float(sys.argv[9]))
         elif sys.argv[1] == 'ri' or sys.argv[1] == 'rg':
             # Parameters order: Alpha,bank,outfile,infile
             fluctuation_run_no_decay(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6],int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11]))
         elif sys.argv[1] == 'bd' or sys.argv[1] == 'co':
             # Parameters order: Alpha,bank,outfile,infile
             fluctuation_run_extinction_DiGraph(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),float(sys.argv[9]))
         elif sys.argv[1] == 'q':
             # activate quarntine at a given time
             fluctuation_run_catastrophe(float(sys.argv[2]), float(sys.argv[3]) ,int(sys.argv[4]), sys.argv[5],
             sys.argv[6],int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11]),
             float(sys.argv[12]),float(sys.argv[13]),sys.argv[14])
