# the discretion of the map in the coordinates
import numpy as np
from coordinationTransformer import CoordinationTransformer as CT
from analytical_geometry import *
import copy
import matplotlib.pyplot as plt
from math import atan2
from collections import deque
import networkx as nx
import heapq
# resolution of the map
R = 1
# the up bound of
EP_ANGLE = 0.01
EP_DISTANCE = 0.5
# infinity
INF = 1e10


def inside(visited, neighbor):
    for i in visited:
        if i == neighbor:
            return True
    return False


def delete_collinear(list_boundary, list_obstacles):
    cur = list_boundary[0]
    n = len(list_boundary)
    list_to_delete = []
    for i in range(1, n):
        if atan2(list_boundary[i].y - cur.y, list_boundary[i].x - cur.x) == \
           atan2(list_boundary[(i+1)%n].y - list_boundary[i].y, list_boundary[(i+1)%n].x - list_boundary[i].x):
            list_to_delete.append(i)
        else:
            cur = list_boundary[i]
    for j in list_to_delete[::-1]:
        list_boundary.pop(j)

    for obstacle in list_obstacles:
        cur = obstacle[0]
        n = len(obstacle)
        list_to_delete = []
        for i in range(1, n):
            if atan2(obstacle[i].y - cur.y, obstacle[i].x - cur.x) == \
                    atan2(obstacle[(i + 1) % n].y - obstacle[i].y,
                          obstacle[(i + 1) % n].x - obstacle[i].x):
                list_to_delete.append(i)
            else:
                cur = obstacle[i]
        for j in list_to_delete[::-1]:
            obstacle.pop(j)


def from_boundary_get_y(list_point, x):
    min_y, max_y = INF, -INF
    for point in list_point:
        if point.x == x:
            min_y = min(min_y, point.y)
            max_y = max(max_y, point.y)
    return min_y, max_y


def check_valid(pb, po, point):
    flag1 = pb.isInPolygon(point)
    flag2 = False
    for obstacle in po:
        if obstacle.isInPolygon(point):
            flag2 = True
            break
    return flag1 and not flag2


def from_obstacles_get_y(obstacles, x):
    min_y, max_y = INF, -INF
    for obstacle in obstacles:
        for point in obstacle:
            if x == point.x:
                min_y = min(min_y, point.y)
                max_y = max(max_y, point.y)
    return min_y, max_y


def calculate_connectivity(slicer):
    connectivity = 0
    last_data = 0
    open_part = False
    connectivity_parts = []
    for i in range(len(slicer)):
        if last_data == 0 and slicer[i] == 1:
            open_part = True
            start_point = i
        elif last_data == 1 and slicer[i] == 0 and open_part:
            open_part = False
            connectivity += 1
            end_point = i
            connectivity_parts.append((start_point, end_point))
        last_data = slicer[i]
    return connectivity, connectivity_parts


def get_adjacent_matrix(parts_left, parts_right):
    """

    :param parts_left:
    :param parts_right:
    :return:
    """
    adjacent_matrix = np.zeros((len(parts_left), len(parts_right)))
    for i in range(len(parts_left)):
        for j in range(len(parts_right)):
            if min(parts_left[i][1], parts_right[j][1]) - max(parts_left[i][0], parts_right[j][0]) > 0:
                adjacent_matrix[i, j] = 1
    return adjacent_matrix


class BoustrophedonCellularDecomposition:
    def __init__(self, filename_boundary, filename_obstacle):
        self.filename = f"{filename_boundary[:-4]}_{filename_obstacle}"
        self.boundary = []
        self.boundary_cnt = 0
        self.obstacles = []
        self.obstacles_cnt = []
        self.aim = []
        self.boundary_simplified = []
        self.obstacles_simplified = []
        self.map = np.empty((0, 0))
        self.map_x = 0
        self.map_y = 0

        with open(f"boundary/{filename_boundary}", 'r') as file:
            lines = file.readlines()
        ref_parts = lines[0].strip().split(',')
        ref_lon, ref_lat = map(float, ref_parts)
        self.transformer = CT(ref_lon, ref_lat)
        for line in lines:
            parts = line.strip().split(',')
            if len(parts) == 2:
                x, y = map(float, parts)
                x, y = self.transformer.geoToCart(x, y)
                self.boundary.append(Point(x, y))
                self.boundary_cnt += 1

        with open(f"obstacle/{filename_obstacle}", 'r') as file2:
            lines2 = file2.readlines()
        obstacle = []
        cnt = 0
        for line in lines2:
            parts = line.strip().split(',')
            if len(parts) == 2:
                x, y = map(float, parts)
                x, y = self.transformer.geoToCart(x, y)
                obstacle.append(Point(x, y))
                cnt += 1
            else:
                self.obstacles.append(obstacle)
                self.obstacles_cnt.append([cnt])
                cnt = 0
                obstacle = []
        self.obstacles.append(obstacle)
        self.obstacles_cnt.append([cnt])

    def graph(self):
        list_boundary_x = []
        list_boundary_y = []
        for i in range(self.boundary_cnt):
            list_boundary_x.append(self.boundary[i].x)
            list_boundary_y.append(self.boundary[i].y)
        list_obstacles_x = []
        list_obstacles_y = []
        for i in range(len(self.obstacles_cnt)):
            list_obstacle_x = []
            list_obstacle_y = []
            for j in range(len(self.obstacles[i])):
                list_obstacle_x.append(self.obstacles[i][j].x)
                list_obstacle_y.append(self.obstacles[i][j].y)
            list_obstacles_x.append(list_obstacle_x)
            list_obstacles_y.append(list_obstacle_y)
        list_aim_x = []
        list_aim_y = []
        for i in range(len(self.aim)):
            list_aim_x.append(self.aim[i].x)
            list_aim_y.append(self.aim[i].y)
        plt.scatter(list_boundary_x, list_boundary_y, label="boundary", c="blue")
        for i in range(len(list_obstacles_x)):
            plt.scatter(list_obstacles_x[i], list_obstacles_y[i], label="obstacles", c="blue")
        plt.scatter(list_aim_x, list_aim_y, c="green")
        boundary_simplified_x = []
        boundary_simplified_y = []
        obstacles_simplified_x = []
        obstacles_simplified_y = []
        for i in range(len(self.boundary_simplified)):
            boundary_simplified_x.append(self.boundary_simplified[i].x)
            boundary_simplified_y.append(self.boundary_simplified[i].y)
        plt.plot(boundary_simplified_x, boundary_simplified_y, c="red")
        for i in range(len(self.obstacles_simplified)):
            obstacle_simplified_x = []
            obstacle_simplified_y = []
            for j in self.obstacles_simplified[i]:
                obstacle_simplified_x.append(j.x)
                obstacle_simplified_y.append(j.y)
            obstacles_simplified_x.append(obstacle_simplified_x)
            obstacles_simplified_y.append(obstacle_simplified_y)
        for i in range(len(obstacles_simplified_x)):
            plt.plot(obstacles_simplified_x[i], obstacles_simplified_y[i], c="red")
        plt.axis("equal")

        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.title('Scatter Plot of boundary and obstacles')

        plt.legend()
        plt.grid(True)
        plt.show()
        plt.imshow(self.map, cmap='viridis')
        plt.colorbar()
        plt.show()

    def deleteCollinear(self):
        list_boundary = self.boundary.copy()
        list_obstacles = copy.deepcopy(self.obstacles)
        current_x = self.boundary[0].x
        current_y = self.boundary[0].y

        current_temp_x = self.boundary[1].x
        current_temp_y = self.boundary[1].y
        current_temp_index = 1
        index_need_to_delete = []
        for i in range(2, self.boundary_cnt):
            if abs(atan2(current_temp_y - current_y, current_temp_x - current_x) -
                   atan2(self.boundary[i].y - current_temp_y, self.boundary[i].x - current_temp_x)) < EP_ANGLE:
                index_need_to_delete.append(current_temp_index)
            else:
                current_x = current_temp_x
                current_y = current_temp_y
            current_temp_y = self.boundary[i].y
            current_temp_x = self.boundary[i].x
            current_temp_index = i
        for i in index_need_to_delete[::-1]:
            list_boundary.pop(i)
        list_need_to_simplified = []
        nums = len(list_boundary)
        for i in range(nums):
            # print(f"{i}:{list_boundary[i].x}, {list_boundary[i].y}")
            if list_boundary[i].distance(list_boundary[(i+1) % nums]) < EP_DISTANCE:
                list_need_to_simplified.append(i)
                # print(f"{i} need to be delete")
        for i in list_need_to_simplified[::-1]:
            list_boundary.pop(i)

        for i in range(len(self.obstacles)):
            n = len(self.obstacles[i])
            # print(n)
            list_need_to_delete = []
            list_need_to_simplify = []
            current_x = self.obstacles[i][0].x
            current_y = self.obstacles[i][0].y
            current_temp_x = self.obstacles[i][1].x
            current_temp_y = self.obstacles[i][1].y
            current_temp_index = 1
            for j in range(2, n):
                if abs(atan2(current_temp_y - current_y, current_temp_x - current_x) -
                       atan2(self.obstacles[i][j].y - current_temp_y, self.obstacles[i][j].x - current_temp_x)) \
                        < EP_ANGLE:
                    list_need_to_delete.append(current_temp_index)
                    # print(current_temp_index)
                else:
                    current_x = current_temp_x
                    current_y = current_temp_y
                current_temp_y = self.obstacles[i][j].y
                current_temp_x = self.obstacles[i][j].x
                current_temp_index = j
            for j in list_need_to_delete[::-1]:
                list_obstacles[i].pop(j)
            nums = len(list_obstacles[i])
            # print(nums)
            for j in range(nums):
                print(f"{j}:{list_obstacles[i][j].x}, {list_obstacles[i][j].y}")
                if list_obstacles[i][j].distance(list_obstacles[i][(j+1) % nums]) < EP_DISTANCE:
                    # print(j)
                    list_need_to_simplify.append((j+1) % nums)
                    # print(f"{i} need to be delete")
            for j in list_need_to_simplify[:nums-1:-1]:
                # print(j)
                list_obstacles[i].pop(j)

        self.obstacles_simplified = list_obstacles
        self.boundary_simplified = list_boundary

    def to_map(self):
        self.pre_process()
        min_x, min_y = float("inf"), float("inf")
        max_x, max_y = float("-inf"), float("-inf")
        for point in self.boundary:
            min_x = min(point.x, min_x)
            min_y = min(point.y, min_y)
            max_x = max(point.x, max_x)
            max_y = max(point.x, max_y)
        self.map_x = min_x
        self.map_y = min_y
        col = int((max_x - min_x)//R+2)
        row = int((max_y - min_y)//R+2)
        self.map = np.ones((row, col))
        # print(f"{col}, {row}")
        list_boundary = []
        list_obstacles = []
        for point in self.boundary:
            p_x, p_y = int((point.x-min_x+R/2)//R), int((point.y-min_y+R/2)//R)
            # print(f"{p_x}, {p_y}")
            if self.map[p_y][p_x] != 0:
                self.map[p_y][p_x] = 0
                list_boundary.append(Point(p_x, p_y))

        for obstacle in self.obstacles:
            list_obstacle = []
            for point in obstacle:
                p_x, p_y = int((point.x - min_x + R/2) // R), int((point.y - min_y + R/2) // R)
                if self.map[p_y][p_x] != 0:
                    self.map[p_y][p_x] = 0
                    list_obstacle.append(Point(p_x, p_y))
            list_obstacles.append(list_obstacle)

        s = Point(18, 1)
        for x in range(2, col):
            min__y, max__y = from_boundary_get_y(list_boundary, x)
            if min__y == 0 and self.map[1][x] != 0:
                s = Point(x, 1)
                break
        # delete_collinear(list_boundary, list_obstacles)
        # self.boundary_simplified = list_boundary
        # self.obstacles_simplified = list_obstacles
        # pb = Polygon(list_boundary)
        # po = []
        # for obstacle in list_obstacles:
        #     poo = Polygon(obstacle)
        #     po.append(poo)
        #
        # while True:
        #     x = random.randint(0, col)
        #     y = random.randint(0, row)
        #     s = Point(x, y)
        #     if check_valid(pb, po, s):
        #         break
        print(s)

        s = (s.x, s.y)

        visited = [s]
        queue = deque()
        queue.append(s)
        print(f"{col}, {row}")
        while queue:
            node = queue.popleft()
            # print(f"||{node}")
            self.map[node[1]][node[0]] = 2
            list_neighbor = []
            p1 = (node[0]+1, node[1])
            if node[0]+1 < col and self.map[node[1]][node[0]+1] != 0:
                list_neighbor.append(p1)
            p2 = (node[0]-1, node[1])
            if node[0]-1 >= 0 and self.map[node[1]][node[0]-1] != 0:
                list_neighbor.append(p2)
            p3 = (node[0], node[1]+1)
            if node[1]+1 < row and self.map[node[1]+1][node[0]] != 0:
                list_neighbor.append(p3)
            p4 = (node[0], node[1]-1)
            if node[1]-1 >= 0 and self.map[node[1]-1][node[0]] != 0:
                list_neighbor.append(p4)
            for neighbor in list_neighbor:
                # print(neighbor)
                if neighbor not in visited:
                    queue.append(neighbor)
                    visited.append(neighbor)
                    # print(f"{neighbor} in visited")
                    # for i in visited:
                    #     print(i, end=" ")
                    # print()

        for x in range(col):
            for y in range(row):
                if self.map[y][x] == 2:
                    self.map[y][x] = 1
                else:
                    self.map[y][x] = 0

    def to_coordinates(self):
        list_x, list_y = [], []
        row, col = self.map.shape
        for x in range(col):
            for y in range(row):
                if self.map[y][x] != 0:
                    list_x.append(self.map_x + x*R)
                    list_y.append(self.map_y + y*R)
        plt.scatter(list_x, list_y, c="red")
        list_boundary_x = []
        list_boundary_y = []
        for i in range(self.boundary_cnt):
            list_boundary_x.append(self.boundary[i].x)
            list_boundary_y.append(self.boundary[i].y)
        list_obstacles_x = []
        list_obstacles_y = []
        for i in range(len(self.obstacles_cnt)):
            list_obstacle_x = []
            list_obstacle_y = []
            for j in range(len(self.obstacles[i])):
                list_obstacle_x.append(self.obstacles[i][j].x)
                list_obstacle_y.append(self.obstacles[i][j].y)
            list_obstacles_x.append(list_obstacle_x)
            list_obstacles_y.append(list_obstacle_y)
        list_aim_x = []
        list_aim_y = []
        for i in range(len(self.aim)):
            list_aim_x.append(self.aim[i].x)
            list_aim_y.append(self.aim[i].y)
        plt.scatter(list_boundary_x, list_boundary_y, label="boundary", c="blue")
        for i in range(len(list_obstacles_x)):
            plt.scatter(list_obstacles_x[i], list_obstacles_y[i], label="obstacles", c="blue")
        plt.scatter(list_aim_x, list_aim_y, c="green")
        boundary_simplified_x = []
        boundary_simplified_y = []
        obstacles_simplified_x = []
        obstacles_simplified_y = []
        for i in range(len(self.boundary_simplified)):
            boundary_simplified_x.append(self.boundary_simplified[i].x)
            boundary_simplified_y.append(self.boundary_simplified[i].y)
        plt.plot(boundary_simplified_x, boundary_simplified_y, c="red")
        for i in range(len(self.obstacles_simplified)):
            obstacle_simplified_x = []
            obstacle_simplified_y = []
            for j in self.obstacles_simplified[i]:
                obstacle_simplified_x.append(j.x)
                obstacle_simplified_y.append(j.y)
            obstacles_simplified_x.append(obstacle_simplified_x)
            obstacles_simplified_y.append(obstacle_simplified_y)
        for i in range(len(obstacles_simplified_x)):
            plt.plot(obstacles_simplified_x[i], obstacles_simplified_y[i], c="red")

        plt.axis("equal")

        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.title('Scatter Plot of boundary and obstacles')

        plt.legend()
        plt.grid(True)
        plt.show()

    # the implementation of the bcd algorithm
    def boustrophedon_cellular_decomposition(self):
        last_connectivity = 0
        last_connectivity_parts = []
        # the color of the cell itself
        current_cell = 1
        current_cells = []
        matrix = copy.deepcopy(self.map)
        row, col = matrix.shape

        for x in range(col):
            current_slice = matrix[:, x]
            connectivity, connectivity_parts = calculate_connectivity(current_slice)
            if last_connectivity == 0:
                current_cells = []
                for i in range(connectivity):
                    current_cells.append(current_cell)
                    current_cell += 1
            elif connectivity == 0:
                current_cells = []
                continue
            else:
                adjacent_matrix = get_adjacent_matrix(last_connectivity_parts, connectivity_parts)
                new_cells = [0]*len(connectivity_parts)
                row1, col1 = adjacent_matrix.shape
                for y in range(row1):
                    if np.sum(adjacent_matrix[y, :]) == 1:
                        # find the first non-0 element's row index in slicer adjacent_matrix[y, :]
                        new_cells[np.argwhere(adjacent_matrix[y, :])[0][0]] = current_cells[y]
                    elif np.sum(adjacent_matrix[y, :]) > 1:
                        for idx in np.argwhere(adjacent_matrix[y, :]):
                            new_cells[idx[0]] = current_cell
                            current_cell += 1

                for x1 in range(col1):
                    if np.sum(adjacent_matrix[:, x1]) > 1:
                        new_cells[x1] = current_cell
                        current_cell = current_cell + 1
                    elif np.sum(adjacent_matrix[:, x1]) == 0:
                        new_cells[x1] = current_cell
                        current_cell = current_cell + 1
                current_cells = new_cells
            for cell, slicer in zip(current_cells, connectivity_parts):
                for y1 in range(slicer[0], slicer[1]):
                    matrix[y1][x] = cell
            last_connectivity = connectivity
            last_connectivity_parts = connectivity_parts
        self.map = matrix
        return matrix, current_cell

    def pre_process(self):
        list_b = []
        n = self.boundary_cnt
        for i in range(self.boundary_cnt):
            list_b.append(self.boundary[i])
            angle = atan2(self.boundary[(i+1) % n][1] - self.boundary[i][1],
                          self.boundary[(i+1) % n][0] - self.boundary[i][0])
            cur_x, cur_y = self.boundary[i][0] + R/2 * cos(angle), self.boundary[i][1] + R/2 * sin(angle)
            while ((cur_x - self.boundary[(i+1) % n][0])**2+(cur_y-self.boundary[(i+1) % n][1])**2)**(1/2) > R/2:
                list_b.append(Point(cur_x, cur_y))
                cur_x, cur_y = cur_x + R/2 * cos(angle), cur_y + R/2 * sin(angle)
        self.boundary = list_b
        list_os = []
        for obstacle in self.obstacles:
            list_o = []
            n = len(obstacle)
            for i in range(n):
                list_o.append(obstacle[i])
                angle = atan2(obstacle[(i+1) % n][1] - obstacle[i][1],
                              obstacle[(i+1) % n][0] - obstacle[i][0])
                cur_x, cur_y = obstacle[i][0] + R/2 * cos(angle), obstacle[i][1] + R/2 * sin(angle)
                while ((cur_x-obstacle[(i+1) % n][0])**2 + (cur_y-obstacle[(i+1) % n][1])**2)**(1/2) > R/2:
                    # print((cur_x, cur_y))
                    list_o.append(Point(cur_x, cur_y))
                    cur_x, cur_y = cur_x + R/2 * cos(angle), cur_y + R/2 * sin(angle)
            list_os.append(list_o)
        self.obstacles = list_os

    def get_path(self, path):
        for i in path:
            self.aim.append(Point(self.map_x+i[0]*R, self.map_y+i[1]*R))
        self.to_coordinates()

    def wirteToFile(self):
        with open(f"routes/{self.filename}", "w") as file:
            for i in range(len(self.aim)):
                temp_longitude, temp_latitude = self.transformer.cartToGeo(self.aim[i].x, self.aim[i].y)
                file.write(f'{temp_longitude},{temp_latitude}\n')


def f(g, cur, end):
    return g + ((cur[0] - end[0])**2 + (cur[1] - end[1])**2)**(1/2)


def check_collision(world_map, x, y):
    return world_map[y][x] == 0


def check(open, neighbor):
    for val, n in open:
        if n == neighbor:
            return True, val
    return False, int(1e10)


def check_type(current, pre):
    if current[0] - pre[0] == 0 and current[1] - pre[1] == 1:
        type1 = 1
    elif current[0] - pre[0] == 1 and current[1] - pre[1] == 1:
        type1 = 2
    elif current[0] - pre[0] == 1 and current[1] - pre[1] == 0:
        type1 = 3
    elif current[0] - pre[0] == 1 and current[1] - pre[1] == -1:
        type1 = 4
    elif current[0] - pre[0] == 0 and current[1] - pre[1] == -1:
        type1 = 5
    elif current[0] - pre[0] == -1 and current[1] - pre[1] == -1:
        type1 = 6
    elif current[0] - pre[0] == -1 and current[1] - pre[1] == 0:
        type1 = 7
    elif current[0] - pre[0] == -1 and current[1] - pre[1] == 1:
        type1 = 8
    else:
        type1 = 9
    return type1


class TSPPython:
    def __init__(self, matrix, cells):
        self.matrix = copy.deepcopy(matrix)
        self.graph = nx.Graph()
        for i in range(1, cells):
            self.graph.add_node(i)
        self.edge_path = {}
        self.cell_path = {}
        self.coverage_path = []
        self.test = []
        self.start_end = {}

    def make_graph(self):
        row, col = self.matrix.shape
        dict_cell = {}
        for key in self.graph.nodes():
            dict_cell[key] = []
        x, y = 0, 0
        upward = True
        while x < col:
            while upward and y < row - 1:
                # self.matrix[y][x] = int(self.matrix[y][x])
                if self.matrix[y][x] == 0:
                    y += 1
                    continue
                else:
                    dict_cell[self.matrix[y][x]].append((x, y))
                    # print(dict_cell[self.matrix[y][x]])
                if x+1 < col and self.matrix[y][x+1] != 0 and self.matrix[y][x] != self.matrix[y][x+1]:
                    self.graph.add_edge(self.matrix[y][x], self.matrix[y][x+1])
                y += 1
            while not upward and y >= 0:
                # self.matrix[y][x] = int(self.matrix[y][x])
                if self.matrix[y][x] == 0:
                    y -= 1
                    continue
                else:
                    dict_cell[self.matrix[y][x]].append((x, y))
                if x+1 < col and self.matrix[y][x+1] != 0 and self.matrix[y][x] != self.matrix[y][x+1]:
                    self.graph.add_edge(self.matrix[y][x], self.matrix[y][x+1])
                y -= 1
            upward = not upward
            x += 1
        start_end = {}
        self.cell_path = dict_cell
        for key in dict_cell:
            value = dict_cell[key]
            start_end[key] = (value[0], value[len(value) - 1])
        self.start_end = start_end
        edges_weight = {}
        edges_path = {}
        for edge in self.graph.edges():
            # print(edge)
            edges_path[edge], edges_weight[edge] = self.A_star(start_end[edge[0]][1], start_end[edge[1]][0])

        nx.set_edge_attributes(self.graph, edges_weight, name='weight')
        self.edge_path = edges_path

    def A_star(self, start, end):
        row, col = self.matrix.shape
        a_path = []
        open_list = [(f(0, start, end), start)]
        # print(open_list[0][1])
        close_list = []
        cost = {start: 0}
        parent = {}
        while open_list[0][1] != end:
            current = heapq.heappop(open_list)
            # print(start, end, current[1])
            close_list.append(current)
            list_neighbor = []
            next_up = (current[1][0], current[1][1] + 1)
            if next_up[1] < col and not check_collision(self.matrix, next_up[0], next_up[1]):
                list_neighbor.append(next_up)
            next_down = (current[1][0], current[1][1] - 1)
            if next_down[1] >= 0 and not check_collision(self.matrix, next_down[0], next_down[1]):
                list_neighbor.append(next_down)
            next_left = (current[1][0] - 1, current[1][1])
            if next_left[0] >= 0 and not check_collision(self.matrix, next_left[0], next_left[1]):
                list_neighbor.append(next_left)
            next_right = (current[1][0] + 1, current[1][1])
            if next_right[0] < row and not check_collision(self.matrix, next_right[0], next_right[1]):
                list_neighbor.append(next_right)
            next_upleft = (current[1][0] - 1, current[1][1] + 1)
            if next_upleft[0] >= 0 and next_upleft[1] < row \
                    and not check_collision(self.matrix, next_upleft[0], next_upleft[1]):
                list_neighbor.append(next_upleft)
            next_upright = (current[1][0] + 1, current[1][1] + 1)
            if next_upright[0] < col and next_upright[1] < row \
                    and not check_collision(self.matrix, next_upright[0], next_upright[1]):
                list_neighbor.append(next_upright)
            next_downleft = (current[1][0] - 1, current[1][1] - 1)
            if next_upleft[0] >= 0 and next_upleft[1] >= 0 \
                    and not check_collision(self.matrix, next_downleft[0], next_downleft[1]):
                list_neighbor.append(next_downleft)
            next_downright = (current[1][0] + 1, current[1][1] - 1)
            if next_upleft[0] < col and next_upleft[1] >= 0 \
                    and not check_collision(self.matrix, next_downright[0], next_downright[1]):
                list_neighbor.append(next_downright)
            for neighbor in list_neighbor:
                type2 = check_type(neighbor, current[1])
                if type2 % 2 == 0:
                    temp_cost = cost[current[1]] + 2**0.5
                else:
                    temp_cost = cost[current[1]] + 1
                is_in_open, cost_n = check(open_list, neighbor)
                is_in_close, cost_v = check(close_list, neighbor)
                if is_in_open and temp_cost < cost[neighbor]:
                    v = (cost_n, neighbor)
                    open_list.remove(v)
                    cost[neighbor] = temp_cost
                    parent[neighbor] = current[1]
                    heapq.heappush(open_list, (f(temp_cost, neighbor, end), neighbor))
                if is_in_close and temp_cost < cost[neighbor]:
                    close_list.remove((cost_v, neighbor))
                    cost[neighbor] = temp_cost
                    parent[neighbor] = current[1]
                    heapq.heappush(open_list, (f(temp_cost, neighbor, end), neighbor))
                if not is_in_open and not is_in_close:
                    cost[neighbor] = temp_cost
                    heapq.heappush(open_list, (f(cost[neighbor], neighbor, end), neighbor))
                    parent[neighbor] = current[1]

        cur = end
        distance = 0
        while cur in parent:
            a_path.insert(0, cur)
            distance += ((cur[0] - parent[cur][0])**2+(cur[1] - parent[cur][1])**2)**(1/2)
            cur = parent[cur]
            if cur == start:
                a_path.insert(0, cur)
                break
        return a_path, distance

    def TSP_solver(self):
        tsp_path = nx.approximation.traveling_salesman_problem(self.graph, cycle=False)
        return tsp_path

    def make_path(self):
        tsp_path = self.TSP_solver()
        print(tsp_path)
        visited = []
        for i in range(len(tsp_path) - 1):
            if tsp_path[i] not in visited:
                visited.append(tsp_path[i])
                self.coverage_path += self.cell_path[tsp_path[i]]
                # if (tsp_path[i], tsp_path[i+1]) in self.edge_path:
                #     self.coverage_path += self.edge_path[(tsp_path[i], tsp_path[i+1])]
                #     self.test.append(self.edge_path[(tsp_path[i], tsp_path[i+1])])
                # else:
                #     self.coverage_path += self.edge_path[(tsp_path[i+1], tsp_path[i])][::-1]
                #     self.test.append(self.edge_path[(tsp_path[i+1], tsp_path[i])][::-1])
                new_path, dist = self.A_star(self.start_end[tsp_path[i]][1], self.start_end[tsp_path[i+1]][0])
                self.coverage_path += new_path
                self.test.append(new_path)
            else:
                new_path, dist = self.A_star(self.start_end[tsp_path[i]][0], self.start_end[tsp_path[i+1]][0])
                self.coverage_path += new_path
                self.test.append(new_path)
                # print(self.edge_path[(tsp_path[i+1], tsp_path[i])][::-1])
        self.coverage_path += self.cell_path[tsp_path[len(tsp_path)-1]]
        return self.coverage_path

    def visualize(self):
        pos = nx.spring_layout(self.graph)
        # edge_weights = nx.get_edge_attributes(self.graph, 'weight')
        nx.draw_networkx_nodes(self.graph, pos, node_color='yellow', node_size=500)
        nx.draw_networkx_edges(self.graph, pos)
        # nx.draw_networkx_edge_labels(self.graph, pos, edge_labels=edge_weights)
        labels = {i: f"Node {i}" for i in self.graph.nodes()}

        nx.draw_networkx_labels(self.graph, pos, labels=labels)
        plt.axis('off')
        plt.show()

    def map(self):
        x, y = zip(*self.coverage_path)
        plt.plot(x, y, color='green')
        plt.xlabel('X')
        plt.ylabel('Y')
        for i in self.test:
            ix, iy = zip(*i)
            plt.plot(ix, iy, color='red')
        for i in self.start_end:
            plt.scatter(self.start_end[i][0][0], self.start_end[i][0][1], color='yellow')
            plt.scatter(self.start_end[i][1][0], self.start_end[i][1][1], color='blue')
        plt.title('Coverage Path Planning')
        plt.show()


if __name__ == "__main__":
    bcd = BoustrophedonCellularDecomposition("not_convex_4_small.csv", "triple_triangle_small.csv")
    bcd.to_map()
    matrix, cells = bcd.boustrophedon_cellular_decomposition()
    print(cells)
    bcd.graph()
    bcd.to_coordinates()
    graph = TSPPython(matrix, cells)
    graph.make_graph()
    graph.visualize()
    coverage_path = graph.make_path()
    graph.map()
    bcd.get_path(coverage_path)
    bcd.wirteToFile()
