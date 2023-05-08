from curses import delay_output
import numpy as np

from functions import get_base_matrix, get_inverse_matrix

def dual_simplex_method(c_vec, a_mtx, b_vec, b_base, iter = 1, a_base_inv = []):
    #1
    if iter == 1:
        for i in range(len(b_base)):
            b_base[i] -= 1
        a_base_inv = np.linalg.inv(get_base_matrix(a_mtx, b_base))

    #print('inverse matrix', a_base_inv)
    #2
    c_base = []
    for i in b_base:
        c_base.append(c_vec[i])
    #print('c_base', c_base)
    #3
    y_vec = (np.matrix(c_base) * a_base_inv).A1
    #print('yt', y_vec)
    #4
    pseudo_plan = a_base_inv.dot(b_vec)
    #print('pseudo', pseudo_plan)
    #5
    k = [0] * len(c_vec)
    j = 0
    for i in b_base:
        k[i] = pseudo_plan[j]
        j += 1

    #print('kt', k)
    
    j_k = -1
    for i in range(len(k)):
        if k[i] < 0:
            j_k = i

    if j_k == -1:
        return k
    
    j_k_ind = b_base.index(j_k)

    #7
    dlt_y_vec = a_base_inv[j_k_ind]

    #//print('delta yt', dlt_y_vec)

    mu = {}
    for i in range(len(c_vec)):
        if i not in b_base:
            mu[i] = np.array(dlt_y_vec).dot(np.array(a_mtx)[:, i])

    #print(mu)
    #8
    excp = True
    for it in mu:
        excp = excp & (mu[it] >= 0)
    #print(excp)
    if excp:
        raise Exception('задача не совместна')
    #9
    sigma = []
    for i in range(len(c_vec)):
        if i not in b_base and mu[i] < 0:
            sigma.append((c_vec[i] - np.array(a_mtx)[:,i].dot(y_vec)) / mu[i])
    
    sigma_0 = min(sigma)
    j_0 = sigma.index(sigma_0)

    b_base[j_k_ind] = j_0

    tmp = get_inverse_matrix(a_base_inv, np.array(np.array(a_mtx)[:, j_0]), j_k_ind)
    
    return dual_simplex_method(c_vec, a_mtx, b_vec, b_base, iter + 1, np.array(tmp))


def lab4():
    c = [-4, -3, -7, 0, 0]
    a = [[-2, -1, -4, 1, 0],
                        [-2, -2, -2, 0, 1]]
    b = [-1, -1.5]
    b_base = [4, 5]

    print(dual_simplex_method(c, a, b, b_base))

if __name__ == '__main__':
    lab4()