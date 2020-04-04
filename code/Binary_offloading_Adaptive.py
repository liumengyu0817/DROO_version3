import scipy.io as sio                     # import scipy.io for .mat file I/
import numpy as np                         # import numpy
from convex import Optimization
from dnn import MemoryDNN
import time
import gc
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"


# 得到所有offloading的可能性
def exhaust_search(K):
    res, path = [], [0] * K

    def Core(index, K):
        if index == K:
            res.append(np.array(path[:]))
            return
        for i in range(2):
            path[index] = i
            Core(index + 1, K)

    Core(0, K)
    return res


def plot_rate( rate_his, rolling_intv = 50, mode='DROO'):
    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    import pandas as pd
    import matplotlib as mpl

    rate_array = np.asarray(rate_his)
    df = pd.DataFrame(rate_his)

    mpl.style.use('seaborn')
    fig, ax = plt.subplots(figsize=(15, 5))
#    rolling_intv = 20

    plt.plot(np.arange(len(rate_array))+1, df.rolling(rolling_intv, min_periods=1).mean(), 'b')
    plt.fill_between(np.arange(len(rate_array))+1, df.rolling(rolling_intv, min_periods=1).min()[0], df.rolling(rolling_intv, min_periods=1).max()[0], color = 'b', alpha = 0.2)
    plt.ylabel('Normalized Computation Rate')
    plt.xlabel('Time Frames')
    # plt.ylabel('归一化代价')
    # plt.xlabel('时间')
    # plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    if mode == 'DROO':
        plt.savefig("../pic/pics/DROO_convergence_ratio_K{0}_D{1}_MonteCarlo{2}.jpg".format(str(N), str(Distance),
                                                                                           str(n)))
    if mode == 'Adaptive':
        plt.savefig("../pic/pics/Adaptive_convergence_ratio_K{0}_D{1}_MonteCarlo{2}.jpg".format(str(N), str(Distance),
                                                                                            str(n)))
    # plt.show()




def save_to_txt(rate_his, file_path):
    with open(file_path, 'w') as f:
        for rate in rate_his:
            f.write("%s \n" % rate)


if __name__ == "__main__":
    '''
        This algorithm generates K modes from DNN, and chooses with largest
        reward. The mode with largest reward is stored in the memory, which is
        further used to train the DNN.
        Adaptive K is implemented. K = max(K, K_his[-memory_size])
    '''
    Distance = 10
    N = 10                     # number of users
    n = 30000                     # number of time frames
    K = N                   # initialize K = N
    decoder_mode = 'OP'    # the quantization mode could be 'OP' (Order-preserving) or 'KNN' or 'exhuast' or 'new_quan'
    Memory = 1024          # capacity of memory structure
    Delta = 10000             # Update interval for adaptive K

    print('#user = %d, #channel=%d, K=%d, decoder = %s, Memory = %d, Delta = %d'%(N,n,K,decoder_mode, Memory, Delta))

    # Load data
    channel = np.loadtxt('../data/K{0}_D{1}/InputData_MonteCarlo{2}_K{0}_Distance{1}.txt'.format(str(N), str(Distance), str(n)))
    optimal_rate = np.loadtxt('../data/K{0}_D{1}/Binary_LR_cost_MonteCarlo{2}_K{0}_Distance{1}.txt'.format(str(N), str(Distance), str(n)))

    # dnn 输入归一化
    channel_h, channel_m, channel_C = channel[:, :2], channel[:, 2:N+2], channel[:, N+2:]
    print('channel_h', channel_h.shape, 'channel_m', channel_m.shape, 'channel_C', channel_C.shape)
    channel_h = (channel_h - channel_h.min()) / (channel_h.max() - channel_h.min())
    channel_m = (channel_m - channel_m.min()) / (channel_m.max() - channel_m.min())
    channel_C = (channel_C - channel_C.min()) / (channel_C.max() - channel_C.min())
    channel_normalize = np.hstack((channel_h, channel_m, channel_C))
    del channel_h, channel_m, channel_C
    gc.collect()

    mem = MemoryDNN(net=[2*N+2, 120, 80, N],
                    learning_rate=0.01,
                    training_interval=16,
                    batch_size=128,
                    memory_size=Memory
                    )
    optimizer = Optimization(0.1)

    rate_his = np.zeros(n)
    rate_his_ratio = np.zeros(n)
    k_idx_his = np.zeros(n)

    print('###############################')
    print(
        '#user = %d, #channel=%d, K=%d, decoder = %s, Memory = %d, Delta = %d' % (N, n, K, decoder_mode, Memory, Delta))
    start_time = time.time()

    for i in range(n):

        if i == 0.8 * n:
            second_time = time.time()

        if i % (n//10) == 0:
           print('episode: ', i/n)

        if i > 0 and i % Delta == 0:
            if sum(k_idx_his[i - Delta:i] == 0) / Delta >= 0.7:
                K = 1
            else:
                K = N



        # 载入当前信道
        h = channel_normalize[i, :]

        # the action selection must be either 'OP' or 'KNN'
        m_list = mem.decode(h, K, decoder_mode)

        r_list = []

        for m in m_list:
            r_list.append(optimizer.optimize(channel[i, :], m)[0])

        rate_his[i] = np.min(r_list)
        rate_his_ratio[i] = rate_his[i]/optimal_rate[i]
        k_idx_his[i] = np.argmin(r_list)

        # encode the mode with largest reward
        mem.encode(h, m_list[np.argmin(r_list)])

    total_time = time.time()-start_time
    later_time = time.time()-second_time
    mem.plot_cost(N, Distance, n, mode='Adaptive')
    plot_rate(rate_his_ratio, mode='Adaptive')

    print('Total time consumed:%s'%(total_time/3600))
    print('Average time per channel:%s' % (total_time / n))
    print('Later time consumed:%s' % later_time)
    print('Later time per channel:%s'%(later_time/0.2/n))

    # save data into txt
    save_to_txt(rate_his,
                "../data/K{1}_D{2}/Binary_Adaptive_cost_MonteCarlo{0}_K{1}_Distance{2}.txt".format(str(n), str(N),
                                                                                                   str(Distance)))

    # save_to_txt(k_idx_his, "../data/K{1}_D{2}/Binary_Adaptive_k_MonteCarlo{0}_K{1}_Distance{2}.txt".format(str(n), str(N), str(Distance)))
    # save_to_txt(rate_his_ratio, "../data/K{1}_D{2}/Binary_Adaptive_CostRatio_MonteCarlo{0}_K{1}_Distance{2}.txt".format(str(n), str(N), str(Distance)))
    # save_to_txt(mem.cost_his, "../data/K{1}_D{2}/Binary_Adaptive_loss_MonteCarlo{0}_K{1}_Distance{2}.txt".format(str(n), str(N), str(Distance)))






