from sklearn.cluster import AgglomerativeClustering
import numpy
import pandas
import time
import itertools
import os
import sys

f = open(sys.argv[1], "r")      # file with list of files to cluster
files = f.readlines()
f.close()

links = ['single', 'average', 'complete']

levels = {'species_cluster' : 0.95, 'genus_cluster': 0.7}

fun = lambda x: os.path.split(x)[1].split('.')[0]

for file in files:
    file = file.rstrip()

    f = open(file, "r")
    
    while True:
        s = f.readline().rstrip()
        if s == "[no_input_sequences]":
            break
        
    s = f.readline().rstrip()
    size = eval(s)

    while True:
        s = f.readline().rstrip()
        if s == "[input_sequences]":
            break

    ANI = numpy.zeros((size, size))

    f_info = []
    columns = {}
    columns["genome"] = []
    
    for i in range(size):
        s = f.readline().rstrip().split(' ')
        columns["genome"].append(s[1].split('.')[0])
        f_info.append(eval(s[2]))

    while True:
        s = f.readline().rstrip()
        if s == "[lz_similarities]":
            break
    
    lines = f.readlines()
    
    for s in lines:
        x = s.rstrip().split(' ')
        
        if len(x) < 6:
            break
        id0 = eval(x[0])
        id1 = eval(x[1])
        
        a = (eval(x[2]) + eval(x[5])) / (f_info[id0] + f_info[id1])
        ANI[id0, id1] = a
        ANI[id1, id0] = a

    numpy.fill_diagonal(ANI, 1)

    print(f"|ANI >= 0.95| = {numpy.sum(ANI >= 0.95)}")
    print(f"|ANI >= 0.70| = {numpy.sum(ANI >= 0.70)}")

    D = 1 - ANI

    for link in links:
        print(f'{link} link', end='')
        t = time.perf_counter()

        for name in levels.keys():
            threshold = 1 - levels[name]
            clustering = AgglomerativeClustering(affinity='precomputed', linkage=link,
                                                 distance_threshold=threshold, n_clusters=None)
            clustering.fit(D)
            columns[name] = clustering.labels_[:]

        dt = time.perf_counter() - t
        print(f' ({dt:.2f} s)')

        out = pandas.DataFrame(data=columns)
        head, tail = os.path.split(file)
        out.to_csv(sys.argv[2] + f'{tail}_{link}.csv', index=False)     # output path

    print('Done!')