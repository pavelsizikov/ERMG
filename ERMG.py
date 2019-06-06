import networkx as nx
import random as rnd


class Graph(nx.Graph):

    def __init__(self):
        nx.Graph.__init__(self)
        self.clusters = {}

    def Print_Clusters(self):
        # вывод кластеров
        print("\nКластеризованная структура графа:")
        for cluster in self.clusters:
            print(str(cluster + 1) + " : " + str(self.clusters[cluster]))

    def Print_Info(self):
        print("\nКоличество вершин в графе :", self.number_of_nodes())
        print("Количество рёбер в графе :", self.number_of_edges())
        print("Количество компонент связности : ", nx.number_connected_components(self))

    def Print_Diam(self):
        diam = self.Diam()
        print("Диаметр: " + str(diam))

    def Print_CliqDistr(self):
        cliq_dist = self.CliqDistr()
        print("\nРаспределение клик: ")
        for cliq in sorted(cliq_dist.keys()):
            print("Размерность: " + str(cliq) + " Количество: " + str(len(cliq_dist[cliq])))

    def Print_DegDistr(self):
        deg_dist = self.DegDistr()
        print("\nРаспределение степеней: ")
        for i in sorted(deg_dist.keys()):
            print("Степень: " + str(i) + " Количество вершин: " + str(deg_dist[i]))

    # вычисление распределения степеней возвращает словарь ключ - степень, значение - число вершин
    def DegDistr(self):
        deg_distr = {}
        for i in self.nodes:
            degree = nx.degree(self, i)
            if degree in deg_distr:
                deg_distr[degree] += 1
            else:
                deg_distr[degree] = 1
        return deg_distr

# вычисление распределения клик
def CliqDistr(g):
    cliq_distr = {2: set(tuple(edge) for edge in g.edges)}
    new_cnt = 2
    while True:
        cnt = new_cnt
        is_added = False
        for cliq in cliq_distr[cnt]:
            is_connected = True
            for node in g.nodes:
                if node not in cliq:
                    for cliq_node in cliq:
                        if not g.has_edge(node, cliq_node):
                            is_connected = False
                            break
                    if is_connected:
                        if not is_added:
                            is_added = True
                            new_cnt += 1
                        if new_cnt not in cliq_distr.keys():
                            cliq_distr[new_cnt] = set()
                        elmt = [node]
                        elmt.extend(cliq)
                        cliq_distr[new_cnt].add(tuple(sorted(elmt)))
        if cnt == new_cnt:
            break
    return cliq_distr


# вычисление диаметра графа
def Diam(g):
    # если нет рёбер, или больше одной компоненты связности, то диаметр = 0
    if g.number_of_edges() == 0 or nx.number_connected_components(g) > 1:
        return 0
    diameter = 0
    for node in g.nodes:
        max_dist = __Calc_Dist(g, node)
        if max_dist > diameter:
            diameter = max_dist
    return diameter


# Вычисление максимального пути для i-ой вершины
def __Calc_Dist(g, i):
    visited = [i]
    current = tuple([i])
    new_depth = 0
    while True:
        depth = new_depth
        new_current = []
        is_added = False
        for node in current:
            neighbours = g.neighbors(node)
            for nbr in neighbours:
                if nbr not in visited:
                    if not is_added:
                        new_depth += 1
                        is_added = True
                    visited.append(nbr)
                    new_current.append(nbr)
        if new_depth == depth:
            break
        current = tuple(new_current)
    return depth


def ERMG(alpha, pi, n):
    # допустимое отклонение для чисел с плавающей точкой(точность)
    eps = 0.0000001
    warning = "Произошла ошибка. Для выхода нажмите enter..."
    # Проверки введённых данных
    for i in range(len(alpha)):
        if alpha[i] > 1 or alpha[i] < 0:
            print("\nalpha_i должны принадлежать промежутку [0..1]")
            input(warning)
            exit(1)

    if abs(sum(alpha) - 1.0) >= eps:
        print("\nСумма вероятностей alpha_i не равна единице")
        input(warning)
        exit(1)

    if len(alpha) != len(pi):
        print("\nРазмерности alpha_i и pi_ij не совпадают")
        input(warning)
        exit(1)
    # Проверка размерности pi_ij
    for i in range(len(pi)):
        for j in range(i, len(pi[i])):
            if len(pi) != len(alpha):
                print("\nРазмерности alpha_i и pi_ij не совпадают")
                input(warning)
                exit(1)
            if pi[i, j] > 1 or pi[i, j] < 0:
                print("\npi_ij должны принадлежать промежутку [0..1]")
                input(warning)
                exit(1)
    # Инициализируем объект графа
    g = Graph()
    # формируем интервалы попадания вероятности для определения кластера
    intervals = []
    curr_sum = 0
    for i in range(len(alpha)):
        curr_sum += alpha[i]
        intervals.append(curr_sum)
        g.clusters[i] = []

        # заносим вершины в кластеры
    for i in range(n):
        p = rnd.random()
        # Реализуем функцию распределения
        for val in intervals:
            if p < val:
                ind = intervals.index(val)
                g.clusters[ind].append(i + 1)
                break

    for i in range(n):
        g.add_node(i + 1)

    # Генерация модели ERMG(рёбер в кластеризованной структуре)
    for cluster in g.clusters:
        for another_cluster in g.clusters:
            if another_cluster >= cluster:
                # для вершин из кластера i
                for i in g.clusters[cluster]:
                    # проверяем ребро к вершине из кластера j (кластеры могут совпадать, но не повторяться)
                    for j in g.clusters[another_cluster]:
                        # исключаем петли и двойную вероятность для рёбер
                        if i != j and (not (cluster == another_cluster and i < j)):
                            p = rnd.random()
                            if p < pi[cluster][another_cluster]:
                                g.add_edge(i, j)
    return g
