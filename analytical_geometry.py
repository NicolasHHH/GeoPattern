# this is the implementation of boustrophedon decomposition
import copy
from math import atan2, cos, sin
import matplotlib.pyplot as plt

# the precision for define the equality of two vectors
EP = 0.000001
EP_PARALLEL = 1e-3
# the parameter of the differential car
L = 2  # [m]
R = L / 2  # radius of the differential car
STEP = 0.1   # [m]


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return abs(self.x - other.x) < EP and (self.y - other.y) < EP

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return f"({self.x}, {self.y})"

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __mul__(self, number):
        return Point(self.x * number, self.y * number)

    def __neg__(self):
        return Point(-self.x, -self.y)

    def __getitem__(self, item):
        if item == 0:
            return self.x
        elif item == 1:
            return self.y
        else:
            raise "ERROR: not x or y coordinate"

    def __setitem__(self, key, value):
        if key == 0:
            self.x = value
        elif key == 1:
            self.y = value
        else:
            raise "ERROR: not x or y coordinate"
        return self.__str__()

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def distance(self, other):
        return ((self.x - other.x) ** 2 + (self.y - other.y) ** 2) ** (1 / 2)

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def cross(self, other):
        return self.x * other.y - self.y * other.x

    def rotate(self, angle):
        new_x = self.x * cos(angle) - self.y * sin(angle)
        new_y = self.x * sin(angle) + self.y * cos(angle)
        self.x = new_x
        self.y = new_y
        return Point(new_x, new_y)

    def foldX(self):
        new_y = -self.y
        self.y = new_y
        return Point(self.x, new_y)

    def foldY(self):
        new_x = -self.x
        self.x = new_x
        return Point(new_x, self.y)


# the angle between the two vector, the former one need to rotate to be the next one
def angleCounterclockwise(vec1, vec2):
    return atan2(vec1.cross(vec2), vec1.cross(vec2))


# the cross product of (point1 - point3) denoted by A and (point2 - point3) denoted by B
# result > 0 A is in the counterclockwise of B
# result == 0 they are collinear
# result < 0 B counterclockwise A
def multiply(point1, point2, point3):
    return (point1[0] - point3[0]) * (point2[1] - point3[1]) - (point1[1] - point3[1]) * (point2[0] - point3[0])


class LineSEG:
    def __init__(self, a, b):
        self.point1 = a
        self.point2 = b

    def online(self, p):
        return ((multiply(self.point1, p, self.point2) == 0) and
                (((p.x - self.point1.x) * (p.x - self.point2.x) <= 0) and
                 ((p.y - self.point1.y) * (p.y - self.point2.y) <= 0)))

    # to find if two line segments are intersected
    def intersect_1(self, other):
        return ((max(self.point1.x, self.point2.x) >= min(other.point1.x, other.point2.x)) and
                (max(other.point1.x, other.point2.x) >= min(self.point1.x, self.point2.x)) and
                (max(self.point1.y, self.point2.y) >= min(other.point1.y, other.point2.y)) and
                (max(other.point1.y, other.point2.y) >= min(self.point1.y, self.point2.y)) and
                (multiply(other.point1, self.point2, self.point1) * multiply(self.point2, other.point2,
                                                                             self.point1) >= 0)
                and
                (multiply(self.point1, other.point2, other.point1) * multiply(other.point2, self.point2,
                                                                              other.point1) >=
                 0))


class Line:
    # Ax+By+C = 0 we assign A > 0
    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C

    def distancePoint2Line(self, point):
        return abs(self.A * point.x + self.B * point.y + self.C) / (self.A ** 2 + self.B ** 2) ** (1 / 2)

    def isParallel(self, other):
        return abs(self.A * other.B - self.B * other.A) < EP_PARALLEL

    def distanceLine2Line(self, other):
        if not self.isParallel(other):
            return False
        if self.A != 0:
            return abs(other.A * self.C - self.A * other.C) / (other.A * (self.A ** 2 + self.B ** 2) ** (1 / 2))
        else:
            return abs(other.B * self.C - self.B * other.C) / (other.B * (self.A ** 2 + self.B ** 2) ** (1 / 2))


def makeLine(point1, point2):
    sign = 1
    A = point2.y - point1.y
    if (A < 0):
        sign = -1
    A = sign * A
    B = sign * (point1.x - point2.x)
    C = sign * (point1.y * point2.x - point1.x * point2.y)
    return Line(A, B, C)


# decide if two line are intersected
def lineIntersect(line1, line2):
    is_intersect = line1.isParallel(line2)
    det = line1.A * line2.B - line1.B * line2.A
    if is_intersect:
        return True, Point((line2.C * line1.B - line1.C * line2.B) / det,
                           (line2.A * line1.C - line1.A * line2.C) / det)
    else:
        return False, Point(0, 0)


# decide if the line segment seg1 intersect with the seg2's line, if so return True
def lineSEGlineSEGlineintersect(seg1, seg2):
    return multiply(seg1.point1, seg2.point2, seg2.point1) * multiply(seg2.point2, seg1.point2, seg2.point1) >= 0


# decide if line segment seg are intersected with line
def lineSEGlineIntersect(seg, line):
    line2 = makeLine(seg.point1, seg.point2)
    is_intersect, point = lineIntersect(line2, line)
    if is_intersect:
        if seg.online(point):
            return True, point
    return False, Point(0, 0)


# build line which is parallel to the y-axis through point
def makeLine2(point):
    A = 1
    B = 0
    C = -point.x
    return Line(A, B, C)


# return True if seg1 intersect with seg2 and the intersect point is not on end point
def intersect_2(seg1, seg2):
    return ((seg1.intersect_1(seg2)) and
            (not seg1.online(seg2.point1)) and
            (not seg1.online(seg2.point2)) and
            (not seg2.online(seg1.point2)) and
            (not seg2.online(seg1.point1)))


class Polygon:
    def __init__(self, point_list):
        # the points should be arranged in counterclockwise order
        # the adjacent points are connected
        self.point_list = copy.deepcopy(point_list)
        self.edge_list = []
        n = len(point_list)
        for i in range(n):
            self.edge_list.append(LineSEG(point_list[i], point_list[(i + 1) % n]))

        self.rotated_angle = 0
        self.folded_x = False
        self.folded_y = False

        self.aim = []
        self.aim_cnt = 0

    def areaOfPolygon(self):
        nums = len(self.point_list)
        if nums < 3:
            return 0
        result = self.point_list[0].y * (self.point_list[nums - 1].x - self.point_list[1].x)
        for i in range(1, nums, 1):
            result += self.point_list[i] * (self.point_list[i - 1].x - self.point_list[(i + 1) % nums].x)
        return result / 2

    def isCounterclockwise(self):
        if self.areaOfPolygon() > 0:
            return True
        return False

    def reverse(self):
        self.point_list = self.point_list[::-1]

    def isSimple(self):
        n = len(self.point_list)
        for i in range(n):
            seg1 = LineSEG(self.point_list[i], self.point_list[(i + 1) % n])
            c = n - 3
            while c:
                seg2 = LineSEG(self.point_list[(i + 2) % n], self.point_list[(i + 3) % n])
                if seg1.intersect_1(seg2):
                    break
                c -= 1
            if c:
                return False
        else:
            return True

    def checkConvex(self):
        index = 0
        tp = self.point_list[0]
        n = len(self.point_list)
        is_convex = [True] * n
        for i in range(1, n, 1):
            if self.point_list[i].y < tp.y or (self.point_list[i].y == tp.y and self.point_list[i].x < tp.x):
                tp = self.point_list[i]
                index = i
        cnt = n - 1
        is_convex[index] = False
        while cnt:
            if multiply(self.point_list[(index + 1) % n], self.point_list[(index + 2) % n],
                        self.point_list[index]) >= 0:
                is_convex[(index + 1) % n] = False
            else:
                is_convex[(index + 1) % n] = True
            index += 1
            cnt -= 1
        return is_convex

    def isConvex(self):
        is_convex = self.checkConvex()
        for i in is_convex:
            if not i:
                return False
        return True

    # return 2 if inside polygon, return 1 in on polygon's edges,  return 0 if outside polygon
    def isInPolygon(self, point):
        seg1 = LineSEG(point, Point(float("inf"), point.y))
        n = len(self.point_list)
        cnt = 0
        for i in range(n):
            seg2 = LineSEG(self.point_list[i], self.point_list[(i + 1) % n])
            if seg2.online(point):
                return 1
            b_intersect = intersect_2(seg1, seg2)
            b_online1 = seg1.online(self.point_list[(i + 1) % n])
            b_online2 = seg1.online(self.point_list[(i + 2) % n])
            r1 = multiply(self.point_list[i], self.point_list[(i + 1) % n], seg1.point1) * \
                 multiply(self.point_list[(i + 1) % n], self.point_list[(i + 2) % n], seg1.point1)
            r2 = multiply(self.point_list[i], self.point_list[(i + 2) % n], seg1.point1) * \
                 multiply(self.point_list[(i + 2) % n], self.point_list[(i + 3) % n], seg1.point1)
            if b_intersect or b_online1 and not b_online2 and r1 > 0 or b_online2 and r2 > 0:
                cnt += 1
        if cnt % 2 == 1:
            return 0
        else:
            return 2

    # return true if inside convex polygon, return false else, must be decide if it's convex polygon
    def isInConvexPolygon(self, point):
        left, right = 1, len(self.point_list) - 2
        while left <= right:
            mid = (left + right) >> 1
            sig1 = multiply(self.point_list[mid], point, self.point_list[0])
            sig2 = multiply(self.point_list[mid + 1], point, self.point_list[0])
            if sig1 >= 0 and sig2 <= 0:
                if multiply(self.point_list[mid + 1], point, self.point_list[mid]) >= 0:
                    return True
                else:
                    return False
            elif sig1 < 0:
                right = mid - 1
            else:
                left = mid + 1
        return False

    def bestAngle(self):
        angle = 0
        width = float("inf")
        n = len(self.point_list)
        for i in range(n):
            t = 1
            line1 = makeLine(self.point_list[i], self.point_list[(i + 1) % n])
            distance = float("-inf")
            for j in range(n):
                if line1.distancePoint2Line(self.point_list[t]) >= distance:
                    distance = line1.distancePoint2Line(self.point_list[t])
                    record = t
                    print(f"distance {distance} : {self.point_list[record]}, {self.point_list[record]}")
                t = (t + 1) % n
            # print(f"({list_x[i]},{list_y[i]}) -> ({list_x[record]},{list_y[record]}): distance:{distance}")
            if width >= distance:
                width = distance
                angle = angleCounterclockwise(Point(0, 1), self.point_list[(i + 1) % n] - self.point_list[i])
            # print(angle)
        return angle, width

    def rotate(self, angle):
        self.rotated_angle = angle
        for i in self.point_list:
            i.rotate(angle)

    def inverseRotate(self):
        for i in self.point_list:
            i.rotate(-self.rotated_angle)
        for j in self.aim:
            j.rotate(-self.rotated_angle)
        self.rotated_angle = 0

    def pathPlanningStartPoint(self):
        best_angle, width = self.bestAngle()
        self.rotate(-best_angle)
        res1 = float("inf")
        res2 = float("-inf")
        for i in self.point_list:
            res1 = min(res1, i.x)
            res2 = max(res2, i.x)
        line1 = makeLine2(res1 + R)
        line2 = makeLine2(res2 - R)
        list_intersect_1 = []
        list_intersect_2 = []
        for i in self.edge_list:
            is_intersect_1, point1 = lineSEGlineIntersect(i, line1)
            is_intersect_2, point2 = lineSEGlineIntersect(i, line2)
            if is_intersect_1:
                list_intersect_1.append(point1)
            if is_intersect_2:
                list_intersect_2.append(point2)
        p1 = Point(res1 + R, min(list_intersect_1))
        p2 = Point(res1 + R, max(list_intersect_1))
        p3 = Point(res2 - R, min(list_intersect_2))
        p4 = Point(res2 - R, max(list_intersect_2))
        self.inversePoint(p1)
        self.inversePoint(p2)
        self.inversePoint(p3)
        self.inversePoint(p4)
        return [p1, p2, p3, p4]

    def foldX(self):
        self.folded_x = True
        self.point_list = self.point_list[::-1]
        for i in self.point_list:
            i.foldX()

    def inverseFoldX(self):
        if self.folded_x:
            for i in self.point_list:
                i.foldX()
            self.point_list = self.point_list[::-1]
            for j in self.aim:
                j.foldX()
            self.aim = self.aim[::-1]

    def foldY(self):
        self.folded_y = True
        self.point_list = self.point_list[::-1]
        for i in self.point_list:
            i.foldY()

    def inverseFoldY(self):
        if self.folded_y:
            for i in self.point_list:
                i.foldY()
            self.point_list = self.point_list[::-1]
            for j in self.aim:
                j.foldY()
            self.aim = self.aim[::-1]

    def inversePoint(self, point):
        point.rotate(-self.rotated_angle)
        if self.folded_x:
            point.foldX()
        if self.folded_y:
            point.foldY()

    def inverseAll(self):
        self.inverseRotate()
        self.inverseFoldX()
        self.inverseFoldY()

    def pathPlanning(self):
        length = 0.0
        start = float("inf")
        end = float("-inf")
        up_bound = float("-inf")
        down_bound = float("inf")
        for i in self.point_list:
            start = min(start, i.x)
            end = max(end, i.x)
            up_bound = max(up_bound, i.y)
            down_bound = min(down_bound, i.y)
        upward = True
        current_x, current_y = start + R, down_bound - 1
        while current_x < end:
            while upward and current_y <= up_bound:
                # print(f"{current_x},{current_y}", end=" ")
                current_y += STEP
                point = Point(current_x, current_y)
                if self.isInPolygon(point):
                    if self.aim_cnt > 1:
                        length += point.distance(self.aim[self.aim_cnt - 1])
                    self.aim.append(point)
                    self.aim_cnt += 1
            while not upward and current_y >= down_bound:
                # print(f"{current_x},{current_y}", end=" ")
                current_y -= STEP
                point = Point(current_x, current_y)
                if self.isInPolygon(point):
                    length += point.distance(self.aim[self.aim_cnt - 1])
                    self.aim.append(point)
                    self.aim_cnt += 1
            # attention::in this case L / STEP = 20
            n = L // STEP  # this is float?
            for i in range(int(n)):
                # print(f"{current_x},{current_y}", end=" ")
                current_x += STEP
                point = Point(current_x, current_y)
                if self.isInPolygon(point):
                    length += point.distance(self.aim[self.aim_cnt - 1])
                    self.aim.append(point)
                    self.aim_cnt += 1
            upward = not upward
        end_point = Point(current_x, current_y)
        self.inverseAll()
        self.inversePoint(end_point)
        return self.aim, length, end_point

    def graphMode(self):
        list_x = []
        list_y = []
        for i in self.point_list:
            list_x.append(i.x)
            list_y.append(i.y)
        plt.plot(list_x, list_y, label="cell", c="purple")
        # plt.axis("equal")
        #
        # plt.xlabel('X Coordinate')
        # plt.ylabel('Y Coordinate')
        # plt.title('Scatter Plot of boundary and obstacles')
        #
        # plt.legend()
        # plt.grid(True)
        # plt.show()


def intersectPolygonLine(polygon, line):
    intersect_point = []
    for edge in polygon.edge_list:
        is_intersect, point = lineSEGlineIntersect(edge, line)
        if is_intersect:
            if point not in intersect_point:
                intersect_point.append(point)
    intersect_point.sort(key=lambda Point: Point.y, reverse=False)
    return intersect_point


def intersect_polygons_line(polygons, point, polygon):
    line = makeLine2(point)
    point_list = []
    for i in polygons:
        for j in i.edge_list:
            is_intersect, point2 = lineSEGlineIntersect(j, line)
            if is_intersect and point2 not in point_list and point2 not in polygon:
                point_list.append(point2)
    point_list.sort(key=lambda Point: abs(Point.y-point.y), reverse=False)
    list_res = point_list[0:2:1].sort(key=lambda Point: Point.y, reverse=False)
    return list_res


