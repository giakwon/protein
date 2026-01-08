import time

DIRS = {
    "U": (0, 1),
    "D": (0, -1),
    "L": (-1, 0),
    "R": (1, 0),
}

def orthogonal_transforms(x, y):
    transforms = set()
    points = [
        (x, y),
        (-y, x),
        (-x, -y),
        (y, -x),
    ]
    for px, py in points:
        transforms.add((px, py))
        transforms.add((-px, py))  # reflect over y-axis
        transforms.add((px, -py))  # reflect over x-axis
        transforms.add((-px, -py)) # reflect over origin
    return transforms

def method_1(n):
    structures = []
    def backtrack(x=0, y=0, structure=[], occupied=set()):
        if len(structure) == n:
            structures.append("".join(structure))
            return
        occupied.add((x, y))
        for direction, (dx, dy) in DIRS.items():
            if (x + dx, y + dy) in occupied:
                continue
            structure.append(direction)
            backtrack(x + dx, y + dy)
            structure.pop()
        occupied.remove((x, y))
    backtrack()
    return structures

def method_2(n):
    if n == 0:
        return []
    structures = []
    def backtrack(x=1, y=0, structure=["R"], occupied=set()):
        if len(structure) == n:
            structures.append("".join(structure))
            return
        points = orthogonal_transforms(x, y)
        occupied.update(points)
        for direction, (dx, dy) in DIRS.items():
            if (x + dx, y + dy) in occupied:
                continue
            structure.append(direction)
            backtrack(x + dx, y + dy)
            structure.pop()
        occupied.difference_update(points)
    backtrack()
    return structures

def bench(func, *args, **kwargs):
    start = time.perf_counter()
    result = func(*args, **kwargs)
    end = time.perf_counter()
    # print(result)
    print(f"{func.__name__} took {end - start} seconds")
    return result


if __name__ == "__main__":
    for i in range(1, 40):
        print(f"For n equal to {i}")
        # bench(method_1, i)
        bench(method_2, i)
    # print(method_1(4))
    print(method_2(4))