import numpy as np
from convex import Optimization
import time
import pandas as pd


start = time.time()

# 用 BCD method 解决 binary offloading
K = 8
MonteCarlo = 10000
Distance = 10

channel = np.loadtxt('../data/K{0}_D{1}/InputData_MonteCarlo{2}_K{0}_Distance{1}.txt'.format(str(K), str(Distance), str(MonteCarlo)))
optimal_rate = np.loadtxt('../data/K{0}_D{1}/Binary_LR_cost_MonteCarlo{2}_K{0}_Distance{1}.txt'.format(str(K), str(Distance), str(MonteCarlo)))

Cost = np.zeros(MonteCarlo)
iteration_number = 0

all_rate = []

for cyc in range(MonteCarlo):

    # 载入信道
    h = channel[cyc, :]

    # 初始化最初的决策
    decision = np.random.randint(0, 2, K)

    count = 0
    rate = []
    while True:
        if count == 0:
            optimizer = Optimization(0.1)
            m = decision.copy()
            rate.append(optimizer.optimize(h, m)[0])
        else:
            rate_tmp = np.zeros(K)
            for idx in range(K):
                optimizer = Optimization(0.1)
                m = decision.copy()
                m[idx] = m[idx] ^ 1
                rate_tmp[idx] = optimizer.optimize(h, m)[0]
            if min(rate_tmp) < rate[-1]:
                rate.append(min(rate_tmp))
                index = np.argmin(rate_tmp)
                decision[index] = decision[index] ^ 1
            else:
                break
        count += 1

    Cost[cyc] = rate[-1]
    iteration_number += (count/MonteCarlo)
    all_rate.append(np.array(rate))

    if cyc % (MonteCarlo//10) == 0:
        print('episode: ', cyc/MonteCarlo, decision, Cost[cyc], optimal_rate[cyc])


np.savetxt('../data/Binary_BCD_cost_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), Cost)

total_time = time.time() - start
print('Total time consumed:%s' % total_time)
print('Average time per channel:%s' % (total_time / MonteCarlo))
print('Total iteration: ', iteration_number)


# all_rate = pd.DataFrame(all_rate)
# all_rate.fillna(all_rate.mean(),inplace=True)
# np.savetxt('../data/Binary_DROO_allCost_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), all_rate.values)





