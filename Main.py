import numpy as np
import ERMG as ermg
import scipy.special as sl
import datetime
import os
import networkx as nx
import random as rnd
import math
import itertools
import csv

# в completed_edges.csv 11940435 рёбер
EDGE_COUNT = 11940435
# процент экспериментов от рассматриваемого числа вершин за одно полное чтение csv-файла (в целях производительности)
EXPERIMENTS_PCT = 0.05


def Build_Subgrph(edges, n):
    src = ermg.Graph()
    # получаем номер строки, для взятия первого ребра(пары вершин), с которых будет построен граф
    num = rnd.randint(0, EDGE_COUNT - 1)
    nbrs = set(edges[num])
    src.add_nodes_from(nbrs)
    nbrs.clear()
    nds = src.nodes
    while src.number_of_nodes() != n:
        for edge in edges:
            # только одна из уже рассмотренных вершин
            if (edge[0] in nds) != (edge[1] in nds):
                nbrs.update([edge[0], edge[1]])
        exp_cnt = math.ceil(len(nbrs) * EXPERIMENTS_PCT)
        for exp in range(exp_cnt):
            i = rnd.randint(0, len(nbrs) - 1)
            node = list(nbrs)[i]
            nbrs.remove(node)
            src.add_node(node)
            if src.number_of_nodes() >= n:
                break
    # добавляем все рёбра для полученного множества вершин
    for edge in edges:
        if (edge[0] in nds) and (edge[1] in nds):
            src.add_edge(edge[0], edge[1])
    return src


# построение модели Эрдёша-Реньи
def Build_ER(src):
    n = src.number_of_nodes()
    er = ermg.Graph()
    max_edges = float(sl.comb(n, 2, exact=True))
    p = float(src.number_of_edges())/max_edges
    nx_er = nx.erdos_renyi_graph(n, p)
    er.add_nodes_from(nx_er.nodes)
    er.add_edges_from(nx_er.edges)
    return er


# построение модели Барабаши-Альберт
def Build_BA(src, n):
    ba = ermg.Graph()
    edges_count = src.number_of_edges()
    m = round(edges_count/n)
    nx_ba = nx.barabasi_albert_graph(n, m)
    ba.add_nodes_from(nx_ba.nodes)
    ba.add_edges_from(nx_ba.edges)
    return ba


# построение модели Боллобаша-Риордана
def Build_BR(src):
    n = src.number_of_nodes()
    er = ermg.Graph()
    return er


# построение модели ERMG
def Build_ERMG(src):
    n = src.number_of_nodes()

    deg_distr = src.DegDistr()

    q = len(deg_distr)

    alpha = [0]*q
    i = 0
    for deg in deg_distr.keys():
        alpha[i] = float(deg_distr[deg])/n
        i += 1

    degrees = list(deg_distr.keys())
    pi = np.zeros(shape=(q, q))
    for i in range(q):
        for j in range(i, q):
            pi[i][j] = getPi(i, j, src, degrees)
            pi[j][i] = pi[i][j]
    ermg_graph = ermg.ERMG(alpha, pi, n)
    return ermg_graph


# получение pi_ij для ERMG
def getPi(i, j, src, degrees):
    p = 0
    # вычисление внутренней вероятности появления ребра
    if i == j:
        nodes = []
        for node in src.nodes:
            if nx.degree(src, node) == degrees[i]:
                nodes.append(node)
        edges_count = 0
        for edge in src.edges:
            if edge[0] in nodes and edge[1] in nodes:
                edges_count += 1
        if len(nodes) > 1:
            p = float(edges_count)/float(sl.comb(len(nodes), 2, exact=True))
    # вычисление вероятности появления ребра между группами
    else:
        nodes_i = []
        nodes_j = []
        for node in src.nodes:
            if nx.degree(src, node) == degrees[i]:
                nodes_i.append(node)
            if nx.degree(src, node) == degrees[j]:
                nodes_j.append(node)
        edges_count_ij = 0
        for edge in src.edges:
            if (edge[0] in nodes_i and edge[1] in nodes_j) or (edge[0] in nodes_j and edge[1] in nodes_i):
                edges_count_ij += 1
        # вероятность = количество рёбер между группами (без рёбер внутри группы) / максимальное количество рёбер в
        # двудольном графе
        if len(nodes_i) > 0 and len(nodes_j) > 0:
            p = float(edges_count_ij) / float(len(nodes_i)*len(nodes_j))
    return p


def Remove_Isolated_Nodes(graph):
    nodes = []
    for node in graph.nodes:
        if nx.degree(graph, node) == 0:
            nodes.append(node)
    graph.remove_nodes_from(nodes)
    return graph


n = int(input("Введите количество вершин в графе n = "))

if n <= 0:
    print("Количество вершин указано некорректно. Должно быть > 0.")
    exit(1)

# читаем файл с рёбрами исходного графа
edges = []
with open('completed_edges.csv', newline='') as csvfile:
    edge_rows = csv.reader(csvfile, delimiter=',')
    for edge in edge_rows:
        edges.append(tuple(edge))

k = 100
# словарь характеристик: ключ - модель, значение - список [Кол-во рёбер, Диаметр, Коэф. класт-ции, Плотность]
av_params = {"SRC": [0.0, 0.0, 0.0, 0.0],
             "ER": [0.0, 0.0, 0.0, 0.0, 0.0],
             "BA": [0.0, 0.0, 0.0, 0.0],
             # "BR": [0.0, 0.0, 0.0, 0.0],
             "ERMG": [0.0, 0.0, 0.0, 0.0, 0.0]}

f = open("parameters.txt", "w+")
i = 1
while i <= k:   # - повторные эксперименты, для вычисления средних характеристик
    # получаем подграф
    src = Build_Subgrph(edges, n)
    av_params["SRC"][0] += nx.number_of_edges(src)
    av_params["SRC"][1] += nx.diameter(src)
    av_params["SRC"][2] += nx.average_clustering(src)
    av_params["SRC"][3] += max(ermg.CliqDistr(src).keys())

    # строим модель Эрдёша-Реньи
    er = max(nx.connected_component_subgraphs(Build_ER(src)), key=len)
    av_params["ER"][0] += nx.number_of_edges(er)
    av_params["ER"][1] += nx.diameter(er)
    av_params["ER"][2] += nx.average_clustering(er)
    av_params["ER"][3] += max(ermg.CliqDistr(er).keys())
    av_params["ER"][4] += nx.number_of_nodes(er)

    # строим модель Барабаши-Альберт
    ba = Build_BA(src, n)
    av_params["BA"][0] += nx.number_of_edges(ba)
    av_params["BA"][1] += nx.diameter(ba)
    av_params["BA"][2] += nx.average_clustering(ba)
    av_params["BA"][3] += max(ermg.CliqDistr(ba).keys())

    # строим модель Боллобаша-Риордана
    '''
    br = Build_BR(src)
    av_params["BR"][0] += nx.number_of_edges(br)
    av_params["BR"][1] += nx.diameter(br)
    av_params["BR"][2] += nx.average_clustering(br)
    av_params["BR"][3] += nx.density(br)'''

    # строим модель ERMG
    ermg_g = max(nx.connected_component_subgraphs(Build_ERMG(src)), key=len)
    av_params["ERMG"][0] += nx.number_of_edges(ermg_g)
    av_params["ERMG"][1] += nx.diameter(ermg_g)
    av_params["ERMG"][2] += nx.average_clustering(ermg_g)
    av_params["ERMG"][3] += max(ermg.CliqDistr(ermg_g).keys())
    av_params["ERMG"][4] += nx.number_of_nodes(ermg_g)

    # записываем параметры
    if i == 10 or i == 25 or i == 50 or i == 75 or i == 100:
        f.write("k = %d\r\n" % i)
        for model in av_params.keys():
            f.write(model + '\r\n')
            av_edges = av_params[model][0]/i
            f.write("Рёбра %f\r\n" % av_edges)
            av_diam = av_params[model][1]/i
            f.write("Диаметр %f\r\n" % av_diam)
            av_coeff = av_params[model][2]/i
            f.write("Коэффициент кластеризации %f\r\n" % av_coeff)
            av_cliq = av_params[model][3]/i
            f.write("Максимальная клика %f\r\n" % av_cliq)
            if model == "ERMG" or model == "ER":
                av_nodes = av_params[model][4] / i
                f.write("Вершин %f\r\n" % av_nodes)
            f.write("---------------------------\r\n")
    i += 1
f.close()
