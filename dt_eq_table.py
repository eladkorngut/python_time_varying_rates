import numpy as np
import bisect
from scipy.integrate import quad
from scipy.optimize import fsolve
import warnings
import sys
import networkx as nx
import netinithomo
import os

warnings.filterwarnings('ignore', 'The iteration is not making good progress')


def temporal_change_function(rate_type):
    if rate_type == 'c':
        Beta = float(np.load('parmeters.npy'))
        # Beta, amplitude, frequency = np.load('parmeters.npy')
        fun = lambda t: Beta
        end_time = 1.0
    elif rate_type == 's':
        Beta, amplitude, frequency = np.load('parmeters.npy')
        fun = lambda t: Beta * (1 + amplitude * np.cos(2 * np.pi * t / frequency))
        end_time = 2*np.pi*frequency
    return fun,end_time


def cycle_time(rate_type):
    if rate_type == 'c':
        return lambda t: t
    elif rate_type == 's':
        return lambda t:t%(2*np.pi)


def dt_eq(fun,Num_inf,Alpha,SI_connections,Total_time,r):
    integrand = lambda t: Num_inf * Alpha + SI_connections * fun(t+Total_time)
    integral_fun_t = lambda tf: quad(lambda t: integrand(t + Total_time), 0, tf)[0]
    fun_rand_time = lambda t: integral_fun_t(t) + np.log(r)
    tau = float(fsolve(fun_rand_time, 1.0))
    return tau


# def intalize_table(fun, N, k0, size_r, num_time, Alpha, time_vec, r_vec):
#     table = np.array([])
#     for si in range(N*k0):
#         np.append(table,np.array([]))
#         for inf in range(N):
#             np.append(table[si],np.array([]))
#             for time in range(num_time):
#                 np.append(table[si][inf],np.array([]))
#                 for r in range(size_r):
#                     np.append(table[si][inf][time], [dt_eq(fun, inf,Alpha, si, time_vec[time], r_vec[r])])
#     return table


def intalize_table(fun, N, k0, size_r, num_time, Alpha, time_vec, r_vec):
    table = np.empty([N*k0,N+1,num_time,size_r])
    for si in range(N*k0):
        for inf in range(N+1):
            for time in range(num_time):
                for r in range(size_r):
                    table[si][inf][time][r]=dt_eq(fun, inf,Alpha, si, time_vec[time], r_vec[r])
    return table


def search_table_for_tau(table,n,si,r,time,time_vec,r_vec):
        t_pos = bisect.bisect_left(time_vec, time)
        r_pos = bisect.bisect_left(r_vec, r)
        return table[si][n][t_pos][r_pos]


def continous_vec(size_t,size_r,end_time):
    time_vec = np.linspace(0.0,end_time,size_t)
    r_vec = np.logspace(-8,0.0,size_r)
    return time_vec,r_vec

def create_table(size_t,size_r,rate_type,Alpha,N,k0):
    fun,end_time = temporal_change_function(rate_type)
    time_vec,r_vec = continous_vec(size_t, size_r, end_time)
    table = intalize_table(fun, N, k0, size_r, size_t, Alpha, time_vec, r_vec)
    with open('table.npy', 'wb') as f:
        np.save(f, table)
    with open('tinx.npy', 'wb') as f:
        np.save(f, time_vec)
    with open('rinx.npy', 'wb') as f:
        np.save(f, r_vec)
    return 0

def actasmain():
    N = 20
    k0 = 4
    size_r = 100
    size_t = 100
    rate_type = 's'
    Alpha = 1.0
    table = create_table(size_t, size_r, rate_type, Alpha,N,k0)
    print('This no love song')


if __name__ == '__main__':
    submit = False
    if submit == True:
        actasmain()
    elif sys.argv[1] == 'thx':
        prog = sys.argv[1]
        size_t = int(sys.argv[2])
        size_r = int(sys.argv[3])
        rate_type = sys.argv[4]
        Alpha = float(sys.argv[5])
        N = int(sys.argv[6])
        k0 = int(sys.argv[7])
        Num_different_networks = int(sys.argv[8])
        graphname = sys.argv[9]
        Lam = float(sys.argv[10])
        dir_path = sys.argv[11]
        parts = int(sys.argv[12])
        bank = int(sys.argv[13])
        Num_inf = int(sys.argv[14])
        Num_inital_conditions = int(sys.argv[15])
        create_table(size_t, size_r, rate_type,Alpha, N, k0)
        for n in range(Num_different_networks):
            # G = nx.random_regular_graph(k, N)
            G = nx.complete_graph(N)
            G = netinithomo.intalize_homo_temporal_graph(G)
            infile = graphname + '_' + str(Lam).replace('.', '') + '_' + str(n) + '.pickle'
            nx.write_gpickle(G, infile)
            outfile = 'o' + str(Lam).replace('.', '')
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py ' + str(prog) + ' ' +
                          str(Alpha) + ' ' + str(bank) + ' ' + str(outfile) + ' ' + str(infile) + ' ' +
                          str(Num_inital_conditions) + ' ' + str(Num_inf) + ' ' + str(n)  + ' ' + str(rate_type))


