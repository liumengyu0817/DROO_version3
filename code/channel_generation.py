import numpy as np

for dd in [100]:

    # 参数设置
    K = 8  # 任务个数
    MonteCarlo = 30000  # 信道个数

    # Richan衰落信道模型
    distance = dd
    path_loss = -3

    print('**** K%d, MonteCarlo%d, distance%d ****' % (K, MonteCarlo, distance))

    h_ray = np.sqrt(2)/2 * np.random.rand(MonteCarlo, 2) + np.sqrt(2)/2 * 1j * np.random.rand(MonteCarlo, 2)
    h_lai = np.sqrt(3/4) + np.sqrt(1/4) * h_ray
    h = 1e-3 * abs(h_lai) * abs(h_lai) * pow(distance, path_loss)

    # 产生task <m,C>
    m = np.random.randint(100, 500, (MonteCarlo, K)) * 1000
    # C = np.random.randint(500, 2000, (MonteCarlo, K))
    C = np.random.randint(1000, 2000, (MonteCarlo, K))

    # 合并
    input_data = np.hstack((h, m, C))
    np.savetxt('../data/K{0}_D{1}/InputData_MonteCarlo{2}_K{0}_Distance{1}.txt'.format(str(K), str(distance),
                                                                                       str(MonteCarlo)), input_data)













