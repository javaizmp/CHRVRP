import numpy
import matplotlib.pyplot as plt
import cplex
from scipy.spatial import distance_matrix
import math
import subprocess
from datetime import datetime
import string
import os.path
import time
from statistics import mean
import random
from random import sample
import time
import os
import matplotlib
import scipy.stats as stats
import numpy as np
import random as rnd

for x in range(0,20):
    ######### BEGIN CREATING INSTANCE ##########
    b1 = 35  # demand seed
    b2 = 12  # supply seed
    b3 = 100 + x  # demand loc seed
    b4 = 14 + x  # supply loc seed
    b5 = 40 + x  # vehicle location seed
    b6 = 40  # additional seed if needed


    # generating Customers
    # [[d,s,v,c]]
    pts = [0,5,10,15,20,25]
    #lst = [[20,2,4,6]]

    experiment = [30,3,8,5]

    d = experiment[0]

    Nd = [i for i in range(0,d)]

    # generating Suppliers

    s = experiment[1]

    Ns = [i for i in range(0, s)]

    # generating Depots
    v = experiment[2]

    Nv = [i for i in range(0, v)]
    capa = int(math.ceil(d/v))

    vehicle_cap = experiment[3]

    inst_type = 'x'
    suploc = 'x'

    setting = 'D'+str(d)+'S'+str(s)+'V'+str(v)+'Y'+str(vehicle_cap)
    initialsol = 0
    # assume two commoditiesS
    # generating demand for each of the customers
    rnd1 = numpy.random
    rnd1.seed(b1)
    #demand1 = [1 for i in Nd]

    # generating supply for each of the suppliers
    rnd2 = numpy.random
    rnd2.seed(b2)
    #s1 = [1000 for i in Ns]

    plt = matplotlib.pyplot

    rnd3 = numpy.random
    rnd3.seed(b3)
    loc_Dx = rnd3.rand(len(Nd)) * 100
    loc_Dy = rnd3.rand(len(Nd)) * 100

    rnd4 = numpy.random
    rnd4.seed(b4)
    iloc_Sx = rnd4.uniform(20,80,10)
    iloc_Sy = rnd4.uniform(20,80,10)
    loc_Sx = iloc_Sx[:s]
    loc_Sy = iloc_Sy[:s]

    rnd5 = numpy.random
    rnd5.seed(b5)
    iloc_Vx = rnd5.uniform(0,100,100)
    iloc_Vy = rnd5.uniform(0,100,100)
    loc_Vx = iloc_Vx[:v]
    loc_Vy = iloc_Vy[:v]
        ## end locations

    T1 = list()
    T2 = list()
    T3 = list()
    master = {}
    for i in Nd:
        T1.append("Customer")
        key = 'D' + str(i)
        master[key] = []
        master[key].append(loc_Dx[i])
        master[key].append(loc_Dy[i])

    for i in Ns:
        T2.append("Supplier")
        key = 'S' + str(i)
        master[key] = []
        master[key].append(loc_Sx[i])
        master[key].append(loc_Sy[i])

    for i in Nv:
        T3.append("Vehicle")
        key = 'V' + str(i)
        master[key] = []
        master[key].append(loc_Vx[i])
        master[key].append(loc_Vy[i])

    TypeList = T1 + T2 + T3
    NodeList = Nd + Ns + Nv
    res = [i + str(j) for i, j in zip(TypeList, NodeList)]

    XList = list(loc_Dx) + list(loc_Sx) + list(loc_Vx)
    YList = list(loc_Dy) + list(loc_Sy) + list(loc_Vy)
    coor = [[i, j] for i, j in zip(XList, YList)]

    ######### END CREATING INSTANCE ##########

    ######### BEGIN CREATING SUB DISTANCE MATRICES ##########
    distds = []
    for i in Nd:
        distds.append([])
        for j in Ns:
            x1 = float(loc_Dx[i])
            x2 = float(loc_Sx[j])
            y1 = float(loc_Dy[i])
            y2 = float(loc_Sy[j])
            d1 = (x1 - x2) * (x1 - x2)
            d2 = (y1 - y2) * (y1 - y2)
            f = numpy.sqrt(d1 + d2)
            distds[i].append(f)
    distvs = []
    for i in Ns:
        distvs.append([])
        for j in Nv:
            x1 = float(loc_Sx[i])
            x2 = float(loc_Vx[j])
            y1 = float(loc_Sy[i])
            y2 = float(loc_Vy[j])
            d1 = (x1 - x2) * (x1 - x2)
            d2 = (y1 - y2) * (y1 - y2)
            f = numpy.sqrt(d1 + d2)
            distvs[i].append(f)

    distdd = []
    for i in Nd:
        distdd.append([])
        for j in Nd:
            x1 = float(loc_Dx[i])
            x2 = float(loc_Dx[j])
            y1 = float(loc_Dy[i])
            y2 = float(loc_Dy[j])
            d1 = (x1 - x2) * (x1 - x2)
            d2 = (y1 - y2) * (y1 - y2)
            f = numpy.sqrt(d1 + d2)
            distdd[i].append(f)
    dens = []
    for i in Nd:
        dens.append(sorted(distdd[i], reverse= False)[1])
    ctres = numpy.mean(dens)
    arr = np.array(distdd)
    ######### END CREATING SUB DISTANCE MATRICES ##########

    def blankplot():
        plt.scatter(loc_Dx[0:], loc_Dy[0:], c='c', s=100)
        plt.scatter(loc_Vx[0:], loc_Vy[0:], c='y', s=100, marker='s')
        plt.scatter(loc_Sx[0:], loc_Sy[0:], c='r', s=100, marker='^')
        """for i in Nd:
            plt.annotate('$%d$' % i, (loc_Dx[i] - 0.5, loc_Dy[i] - 1))
        for i in Ns:
            plt.annotate('$%d$' % i, (loc_Sx[i] - 0.5, loc_Sy[i] - 1))
        for i in Nv:
            plt.annotate('$%d$' % i, (loc_Vx[i] - 0.5, loc_Vy[i] - 1))"""
        plt.axis('scaled')
        plt.xlim([0, 100])
        plt.ylim([0,100])
        plt.show()

    #blankplot()

    ######## BEGIN CREATING MASTER INSTANCE FILE WITH TIMESTAMP ##########
    stamp = str(datetime.now().timestamp())
    stamp = stamp.replace('.','')

    def writeinstance():
        mastername = ('instance' + '_size_' + str(d) + '_inst_' + str(x) + '.tsp')
        file = open(mastername, "w")
        file.write('NAME : CHRVRP\n')
        file.write('COMMENT : master instance for CHRVRP\n')
        file.write('EDGE_WEIGHT_TYPE : EUC_2D\n')
        file.write('NODE_COORD_SECTION\n')
        for k in Nv:
            file.write('V' + str(k + 100) + ' ')
            file.write(format(loc_Vx[k],'.3f') + 'e+03 ')
            file.write(format(loc_Vy[k],'.3f') + 'e+03')
            file.write('\n')

        for k in Ns:
            file.write('S' + str(k) + ' ')
            file.write(format(loc_Sx[k],'.3f') + 'e+03 ')
            file.write(format(loc_Sy[k],'.3f') + 'e+03')
            file.write('\n')

        for k in Nd:
            file.write('D' + str(k + 10) + ' ')
            file.write(format(loc_Dx[k],'.3f') + 'e+03 ')
            file.write(format(loc_Dy[k],'.3f') + 'e+03')
            file.write('\n')

        file.write('EOF\n')
        file.close()
    writeinstance()

######## END CREATING MASTER INSTANCE FILE WITH TIMESTAMP ##########