from time import time

for n in (int(1e4), int(1e5), int(1e6)):
    print('---> {:.2E} <---'.format(n))
    d = dict()
    l = list()
    k = range(n)
    for i in k:
        d[i] = None
        l.append(i)

    start = time()
    for i in k:
        _ = i in d
    end = time()
    time_d = (end - start)
    print('dict time: {:.2E}'.format(time_d))
    # start = time()
    # for i in k:
    #     _ = i in l
    # end = time()
    # time_l = (end - start)
    # print('list time: {:.2E}'.format(time_l))




