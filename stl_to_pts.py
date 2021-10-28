#!/usr/bin/env python
"""
More complex script that converts locations in 3D space into a Blueprint
dense filled matrix/block types.
For now, only a Steel cube block is used.

Also implements reading from STL files, and performing the point refinement.
"""

import os
import sys
import math
import time
import struct
import zipfile
import multiprocessing
from copy import copy

# Maximum number of points to attempt to generate per process, for memory bounding
# purposes.
MAX_POINTS_PER_PROCESS = 2000.0

# Build the list of unit vetors
UNIT_VECTORS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1),
                (0, 0, -1)]
# Valid slopes, expressed as 1/m = the number of blocks required to complete the slope.
VALID_SLOPES = [1, 2]

TUPLE_DOT = lambda t1, t2: sum([a * b for a, b in zip(t1, t2)])
TUPLE_ADD = lambda t1, t2: (t1[0] + t2[0], t1[1] + t2[1], t1[2] + t2[2])
TUPLE_SUB = lambda t1, t2: (t1[0] - t2[0], t1[1] - t2[1], t1[2] - t2[2])
TUPLE_MUL = lambda t1, t2: (t1[0] * t2[0], t1[1] * t2[1], t1[2] * t2[2])
TUPLE_LE = lambda t1, t2: (t1[0] < t2[0]) and (t1[1] < t2[1]) and (t1[2] < t2[
    2])
TI = lambda a, t: a[t[0]][t[1]][t[2]]
TUPLE_SCALE = lambda a, t: (a * t[0], a * t[1], a * t[2])

SIGN_S = lambda s: -1 if s < 0 else 1 if s > 0 else 0
SIGN_V = lambda v: tuple([SIGN_S(c) for c in v])


def leq(a, b):
    """
    Given two values, return a boolean or None, depending on whether a < b,
    with a == b returning None.
    """
    if a < b:
        return True
    elif a > b:
        return False
    else:
        return None


class Triple(object):
    """
    Represents a triple of any object type. Is used to represent a point as well
    as a triangle of points.
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __getitem__(self, i):
        if i == 0:
            return self.x
        if i == 1:
            return self.y
        if i == 2:
            return self.z

    def __str__(self):
        return "<%s,%s,%s>" % (str(self.x), str(self.y), str(self.z))

    def __repr__(self):
        return str(self)

    def hexsect(self):
        """
        Given a triangle, partition it into six new triangles using the midpoints
        of each edge and the centroid as new vertices.
        """
        centroid = vmean3(self.x, self.y, self.z)
        mp1 = vmean2(self.x, self.y)
        mp2 = vmean2(self.x, self.z)
        mp3 = vmean2(self.y, self.z)
        return [
            Triple(self.x, mp1, centroid),
            Triple(mp1, self.y, centroid),
            Triple(self.x, mp2, centroid),
            Triple(mp2, self.z, centroid),
            Triple(self.y, mp3, centroid),
            Triple(mp3, self.z, centroid),
        ]

    def shift(self, offset):
        self.x = Triple(*TUPLE_ADD(self.x, offset))
        self.y = Triple(*TUPLE_ADD(self.y, offset))
        self.z = Triple(*TUPLE_ADD(self.z, offset))

    def reflect(self, dim):
        scalar = [1 if i != dim else -1 for i in range(1, 3 + 1)]
        if isinstance(self.x, Triple):
            return Triple(Triple(*TUPLE_MUL(self.x, scalar)),
                          Triple(*TUPLE_MUL(self.y, scalar)),
                          Triple(*TUPLE_MUL(self.z, scalar)))
        else:
            return Triple(*TUPLE_MUL((self.x, self.y, self.x), scalar))


class STLFile(object):
    """
    Represents basic knowledge of the STL file format.
    """
    @staticmethod
    def read_binary_stl(file_descriptor):
        """
        Read all triangles from a binary STL file.
        """
        _ = file_descriptor.read(80)  # Read the global file header
        ntris = struct.unpack("<L", file_descriptor.read(4))[0]

        tris = []
        for _ in range(ntris):
            _ = struct.unpack("<fff", file_descriptor.read(4 * 3))  # Normal
            v1 = struct.unpack("<fff", file_descriptor.read(4 * 3))
            v2 = struct.unpack("<fff", file_descriptor.read(4 * 3))
            v3 = struct.unpack("<fff", file_descriptor.read(4 * 3))
            _ = struct.unpack("<H", file_descriptor.read(2))  # Attributes
            tris.append(Triple(Triple(*v1), Triple(*v2), Triple(*v3)))
        return tris

    @staticmethod
    def read_ascii_stl(file_descriptor):
        """
        Reads a single solid section from an ASCII STL file, returning a list
        marging the triangles from all solids found in the file.
        """
        _solids = []
        solid_name, solid_triangles = STLFile.solid_to_triangles(
            file_descriptor)
        while solid_name is not None:
            sys.stderr.write("Reading solid: %s\n" % solid_name)
            _solids.append((solid_name, solid_triangles))
            solid_name, solid_triangles = STLFile.solid_to_triangles(
                file_descriptor)

        triangles = [e for s in _solids for e in s[1]]
        return triangles

    @staticmethod
    def read_triangles(file_descriptor):
        """
        Given a STL file, determine if it is ASCII or binary, and return a list
        of triangles representing all triangles in the file. For ACSII files,
        the triangles of all solids in the file.
        """
        # First, determine the file is binary or ASCII. Do this by reading the
        # first 5 bytes and if those are 'solid', then assume it is ASCII. Seek
        # back to the start after reading the header.
        is_ascii = (file_descriptor.read(5) == 'solid')
        file_descriptor.seek(0)

        if is_ascii:
            return STLFile.read_ascii_stl(file_descriptor)
        else:
            return STLFile.read_binary_stl(file_descriptor)

    @staticmethod
    def solid_to_triangles(file_descriptor):
        """
        Read through a Solid section of an STL file and extract all triangles as
        Triples of Triples
        """
        # Read in the 'solid' line, getting the name (which may be an empty string)
        # If there's a null in the name, then this is a binary file, and it should
        # be discarded and ignored.
        name_line = file_descriptor.readline()
        if name_line == "" or '\0' in name_line:
            return (None, None)

        name = name_line.strip().split(" ", 1)[1]

        triangles = []
        # Can't mix iteration and .readline(), which made conditionally reading
        # multiple lines awkward, so the infinite while loop and break is used.
        while True:
            line = file_descriptor.readline()

            # If the end of the solid is reached, leave.
            if line.strip().startswith("endsolid"):
                break
            # Otherwise, if it's a vertex line (skipping the loop/endloop lines
            # completely)
            elif line.strip().startswith("vertex"):
                triangles.append(
                    Triple(
                        Triple(
                            *[float(i) for i in line.strip().split(" ")[1:]]),
                        Triple(*[
                            float(i) for i in
                            file_descriptor.readline().strip().split(" ")[1:]
                        ]),
                        Triple(*[
                            float(i) for i in
                            file_descriptor.readline().strip().split(" ")[1:]
                        ])))
        return (name, triangles)


def triangle_list_bounds(tris):
    """
    Given a list of triangles, find the minimum and maximum bounds in each
    dimension.
    """
    bounds = []
    for i in range(3):
        coords = [a[b][i] for a in tris for b in range(3)]
        bounds.append((min(coords), max(coords)))

    return bounds


def vsub(u, v):
    """
    Vector subtraction.
    """
    return Triple(u.x - v.x, u.y - v.y, u.z - v.z)
    #return tuple([i-j for i,j in zip(u,v)])


def l2_norm(v):
    """
    Return the L2 norm of the vector.
    """
    return math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
    #return math.sqrt(sum([i*i for i in v]))


def max_edge_norm(Tri):
    """
    Given an triangle, find length of the longest edge.
    """
    return max(l2_norm(vsub(Tri.x, Tri.y)), l2_norm(vsub(Tri.y, Tri.z)),
               l2_norm(vsub(Tri.x, Tri.z)))


def mean(v):
    """
    Find the mean of all elements of a triple.
    """
    if len(v) > 0:
        return (v.x + v.y + v.z) / 3.0
    else:
        raise ValueError('Cannot compute mean of zero-length vector.')


def vmean2(P1, P2):
    return Triple((P1.x + P2.x) / 2.0, (P1.y + P2.y) / 2.0,
                  (P1.z + P2.z) / 2.0)


def vmean3(P1, P2, P3):
    return Triple((P1.x + P2.x + P3.x) / 3.0, (P1.y + P2.y + P3.y) / 3.0,
                  (P1.z + P2.z + P3.z) / 3.0)


def split_tri(Tri, Resolution):
    """
    Given a single triangle and a spatial resolution, iteratively hexsect the triangles
    so that no subtriangle has an edge length longer than the given resolution.

    Recursive might have been nicer to look at, but this has efficiency benefits for
    large triangles
    """
    small = []
    large = [Tri]
    while len(large) > 0:
        tris = []
        for t in large:
            tris.extend(t.hexsect())
        large = []
        for t in tris:
            if max_edge_norm(t) > Resolution:
                large.append(t)
            else:
                small.append(t)
    return small


def rescale_round_point(Point, Resolution):
    """
    Transform a point to the nearest lattice point that lies on a grid with
    the specified resolution.

    Returned as a tuple for use in a set() for reduction to unique points only.
    """
    return (int(round(Point.x / Resolution)), int(round(Point.y / Resolution)),
            int(round(Point.z / Resolution)))


def parallel_split_tris(Primitives, Resolution, BatchSize=100):
    """
    Perform the split_tris() operation on chunks of primitives in parallel, and
    recombine at the end.
    """
    # For the number of jobs per process, look the bounds, the primitive count,
    # and the resolution. Create process that will have an approximate bound
    # on the number of points generated.
    size = [i[1] - i[0] for i in triangle_list_bounds(Primitives)]
    # Primitives per unit cube of volume, approximately.
    prims_per_unit3 = len(Primitives) / (size[0] * size[1] * size[2])
    # Resolution is essentially units-per-point, so dimensional analysis gives...
    points_per_prim = 1.0 / math.pow(
        (Resolution**3 * prims_per_unit3), 1.0 / 3)
    prims_per_process = int(math.ceil(MAX_POINTS_PER_PROCESS /
                                      points_per_prim))

    # The prims_per_process should be at most enough so that there are more
    # processes than CPUs.
    prims_per_process = min(
        int(math.floor(len(Primitives) / (3 * multiprocessing.cpu_count()))),
        prims_per_process)
    sys.stderr.write(
        "Approximate number of points generated per process: %s\n" %
        (points_per_prim * prims_per_process))

    primitive_chunks = [
        Primitives[i:i + prims_per_process]
        for i in range(0, len(Primitives), prims_per_process)
    ]
    output_queue = multiprocessing.Queue()
    procs = [
        multiprocessing.Process(target=split_tris,
                                args=(chunk, Resolution, BatchSize,
                                      output_queue))
        for chunk in primitive_chunks
    ]
    sys.stderr.write("Prepared %d processes of work\n" % len(procs))

    # First, start cpu_count() processes.
    running_procs = procs[:multiprocessing.cpu_count()]
    for p in running_procs:
        p.start()
    queued_procs = [p for p in procs if p not in running_procs]
    finished_procs = 0

    pts = set()
    # As long as there's a running process, keep cycling.
    while len(running_procs) > 0:
        # Attempt to join all running processes.
        for p in running_procs:
            p.join(0.0)

        while not output_queue.empty():
            pipe_pts = output_queue.get()
            finished_procs += 1
            sys.stderr.write("%d (%d/%d) " %
                             (len(pipe_pts), finished_procs, len(procs)))
            pts.update(pipe_pts)

        # Rebuild the running processes list to only include those still alive
        running_procs = [p for p in running_procs if p.is_alive()]

        # If there are fewer running processes than available CPUs, start some more.
        pending_procs = queued_procs[:multiprocessing.cpu_count() -
                                     len(running_procs)]
        for p in pending_procs:
            p.start()
        # Once started, add them to the pending procs, and remove them from the
        # list of queued processes.
        running_procs += pending_procs
        queued_procs = [p for p in queued_procs if p not in running_procs]

        # Give the processes another second to do some work before checking on them.
        # Prevents a certain amount of busywaiting.
        time.sleep(1.0)

    while not output_queue.empty():
        pts.update(output_queue.get())
        finished_procs += 1
        sys.stderr.write("%d (%d/%d) " %
                         (len(pipe_pts), finished_procs, len(procs)))
    sys.stderr.write("\n")

    return list(pts)


def split_tris(Primitives, Resolution, BatchSize=100, OutputQueue=None):
    """
    Given a list of triangles, split all triangles to the given resolution, and
    flatten the resulting list of triangles to a list of points with duplicates removed.

    Consider batches of triangles together, aggregating the batches of partitioned triangles
    into the cumulative points list, removing duplicates with progression through the
    triangle list. This limits the growth rate of the point and triangle list, significantly
    improving memory consumption.
    """
    start_time = time.time()
    last_print_time = time.time()
    pts = set()
    tris = []
    tris_handled = 0
    for p in Primitives:
        tris.extend(split_tri(p, Resolution))
        # For memory efficiency, perform batch-frequency flatten/union operations, to keep
        # the list of points at any given point in time bounded and reasonable.
        tris_handled += 1
        if (tris_handled % BatchSize) == 0:
            if OutputQueue is None and time.time() - last_print_time > 0.5:
                last_print_time = time.time()
                sys.stderr.write("%d/%d (ETA: %f)\n" %
                                 (tris_handled, len(Primitives),
                                  (len(Primitives) - tris_handled) *
                                  (time.time() - start_time) / tris_handled))
            # sys.stderr.write("Batch done (%d)\n" % tris_handled)
            pts.update([
                rescale_round_point(t[i], Resolution) for i in range(3)
                for t in tris
            ])
            tris = []

    # sys.stderr.write("Final round (%d)\n" % tris_handled)
    # One final round of flatten/union
    pts.update([
        rescale_round_point(t[i], Resolution) for i in range(3) for t in tris
    ])
    pts_l = list(pts)

    # LEGACY: Super slow on pypy (2x CPython), included for posterity and entertainment.
    #pts = list(set([ rescale_round_point(t[i], Resolution) for i in range(3) for t in tris ]))
    if OutputQueue is not None:
        OutputQueue.put(pts_l)
    return pts_l


def p_norm(coords, p):
    """
    Return the p-norm of the list.
    """
    return math.pow(sum([math.pow(c, p) for c in coords]), 1.0 / p)


def lattice_ball(radius, norm=lambda x, y, z: p_norm((x, y, z), 2)):
    """
    Given a radius, find all integer coordinates within that radius of the origin
    in three dimensional Euclidean space (by default).
    """
    coord_range = range(-radius, radius + 1)
    brush = [(x, y, z) for x in coord_range for y in coord_range
             for z in coord_range if norm(x, y, z) <= radius]
    return brush


def parallel():
    shm_stat = None
    try:
        shm_stat = os.stat('/dev/shm')
    except OSError as e:
        if e.errno != 2:
            raise e
    return shm_stat is not None


def parallel_hollow(pts, radius=1):
    """
    Perform model hollowing in parallel across cpu_count() processes.
    """
    return list_parallelize(pts, (radius, pts), hollow)


def hollow(pts, radius=1, all_pts=None, output_queue=None):
    """
    Perform inverted morphological erosion, by keeping all blocks that would have
    failed the erosion test, and discarding all blocks that would have passed the
    erosion test.
    """
    if all_pts is None:
        all_pts = pts

    brush = lattice_ball(radius)
    if len(brush) == 1:
        return pts

    start_time = time.time()
    last_print_time = time.time()
    npts = 0
    orig_pts = dict([(p, False) for p in all_pts])
    for p in orig_pts.keys():
        for b in brush:
            t = TUPLE_ADD(p, b)
            if t not in orig_pts:
                orig_pts[p] = True
                break
        npts += 1
        if output_queue is None and time.time() - last_print_time > 0.5:
            last_print_time = time.time()
            sys.stderr.write("%d/%d (ETA: %f)\n" %
                             (npts, len(pts), (len(pts) - npts) *
                              (time.time() - start_time) / npts))

    ret = [p for p, c in orig_pts.items() if c]
    if output_queue is not None:
        output_queue.put(ret)
    return ret


def list_parallelize(items, args, func):
    """
    Given a list of items, some additional arguments, and a function to call,
    call that function like map but with chosen arguments using Process()
    objects. Assume that the function given takes in a Queue() as a final argument.
    """
    items_per_proc = int(
        math.ceil(1.0 * len(items) / multiprocessing.cpu_count()))
    item_chunks = [
        items[i:i + items_per_proc]
        for i in range(0, len(items), items_per_proc)
    ]
    output_queue = multiprocessing.Queue()
    procs = [
        multiprocessing.Process(target=func,
                                args=(chunk, ) + args + (output_queue, ))
        for chunk in item_chunks
    ]

    for p in procs:
        p.start()
    running_procs = [p for p in procs]

    ret = set()
    while len(running_procs) > 0:
        for p in running_procs:
            p.join(0.0)
        while not output_queue.empty():
            ret.update(output_queue.get())
        time.sleep(0.25)
        running_procs = [p for p in running_procs if p.is_alive()]

    while not output_queue.empty():
        ret.update(output_queue.get())

    return list(ret)


def parallel_morphological_dilate(pts, radius=2):
    """
    Perform morphological dilation in parallel across cpu_count() processes.
    """
    return list_parallelize(pts, (radius, ), morphological_dilate)


def parallel_morphological_erode(pts, radius=2):
    """
    Perform morphological erosion in parallel across cpu_count() processes.
    """
    return list_parallelize(pts, (radius, pts), morphological_erode)


def morphological_dilate(pts, radius=2, output_queue=None):
    """
    Given a list of tuples of integer coordinates, all of the same dimension,
    dilate the list of points to include all points within the given radius of
    a point in the input list.
    """
    brush = lattice_ball(radius)
    if len(brush) == 1:
        return pts
    new_pts = set()

    start_time = time.time()
    last_print_time = time.time()
    npts = 0
    for p in pts:
        for b in brush:
            t = TUPLE_ADD(p, b)
            new_pts.update([t])
        npts += 1
        if output_queue is None and time.time() - last_print_time > 0.5:
            last_print_time = time.time()
            sys.stderr.write("%d/%d (ETA: %f)\n" %
                             (npts, len(pts), (len(pts) - npts) *
                              (time.time() - start_time) / npts))

    new_pts.update(pts)
    ret = list(new_pts)
    if output_queue is not None:
        output_queue.put(ret)
    return ret


def morphological_erode(pts, radius=2, all_pts=None, output_queue=None):
    """
    Given a list of tyuples of integer coordinates, all of the same dimension,
    erode the list of points to include only those points that include every
    point within radius units of it.
    """
    if all_pts is None:
        all_pts = pts

    brush = lattice_ball(radius)
    if len(brush) == 1:
        return pts

    start_time = time.time()
    last_print_time = time.time()
    npts = 0
    orig_pts = dict([(p, True) for p in all_pts])
    for p in pts:
        for b in brush:
            t = TUPLE_ADD(p, b)
            if t not in orig_pts:
                orig_pts[p] = False
                break
        npts += 1
        if output_queue is None and time.time() - last_print_time > 0.5:
            last_print_time = time.time()
            sys.stderr.write("%d/%d (ETA: %f)\n" %
                             (npts, len(pts), (len(pts) - npts) *
                              (time.time() - start_time) / npts))

    ret = [p for p in pts if orig_pts[p]]
    if output_queue is not None:
        output_queue.put(ret)
    return ret


def adjacency_vectors(position, forward, points):
    """
    Return the vector that points out of what should be the bottom of any
    sloped blocks placed.
    """
    # To start, for each unit vector that is perpendicular to the forward vector,
    # check to see if there is a block 'next' to the path. If there is precisely
    # one vector for which this is true, take the director that dots to -1 with
    # the vector satisfying this criteria
    adj = []
    for v in UNIT_VECTORS:
        if TUPLE_DOT(forward, v) != 0:
            continue
        else:
            p = TUPLE_ADD(position, v)
            if p in points and points[p] == 0:
                adj.append(v)
    return [(vec, [v for v in UNIT_VECTORS if TUPLE_DOT(v, vec) == -1][0])
            for vec in adj]


def blocks_to_csv(blocks):
    return [
        "{:d},{:d},{:d},{:d},{:d}".format(b[0], b[1], b[2], b[3], b[4])
        for b in blocks
    ]


def dbm_bitmask(dbm, block_type="\x87"):
    block_strings = ""
    bm = []

    # Current mask, and the bit being twiddled.
    # When this rolls over, it gets appended to the bm list.
    cm = 0
    cb = 1

    for i in range(len(dbm)):
        x = dbm[i]
        for j in range(len(x)):
            y = x[j]
            for k in range(len(y)):
                z = y[k]
                if z != False:
                    #block_strings += "\x87\x01\x00\x00"
                    block_strings += block_type + \
                                     (chr(z[1]) if len(z) > 1 else chr(1)) + \
                                     chr(0) + \
                                     (chr(z[0]) if len(z) > 0 else chr(0))
                    cm += cb
                cb *= 2
                if cb > 128:
                    bm.append(cm)
                    cb = 1
                    cm = 0
    if cb > 1:
        bm.append(cm)
    return (bm, block_strings)


def sparse_to_dense(positions, l, w, h):
    # Map the numeric axes, in terms of the major-ordering of the arrays, to the named axes
    Z = (l, 0)
    Y = (w, 1)
    X = (h, 2)

    a = [[[False for _ in range(Z[0])] for _ in range(Y[0])]
         for _ in range(X[0])]

    for p in positions:
        P = [p[X[1]], p[Y[1]], p[Z[1]]]

    return a


def list_subtract(l1, l2):
    return [l1[i] - l2[i] for i in range(len(l1))]


def bounding_box(positions):
    m = [min([p[i] for p in positions]) for i in range(3)]
    M = [max([p[i] for p in positions]) for i in range(3)]
    return (m, M)


def iterative_flood_fill(dbm, start, VisitedType):
    """
    Use a manual stack context to iteratively flood-fill the volume around the hull.
    """
    M = (len(dbm), len(dbm[0]), len(dbm[0][0]))
    directions = UNIT_VECTORS

    cur_pos = start
    # Each item in the rail is a position, and the directions we have left to try starting at that
    # position.
    trail = [(cur_pos, copy(UNIT_VECTORS))]

    # Prime the pump by testing the start, which is guaranteed to be in-bounds, but not guaranteed
    # to be unoccupied or unvisited. If the position has been visited, but is unoccupied, then the
    # value will be changed from False to None. If the position is occupied, then it will have a
    # value that is neither False nor None.
    if dbm[cur_pos[0]][cur_pos[1]][cur_pos[2]] == None:
        return 0
    elif isinstance(dbm[cur_pos[0]][cur_pos[1]][cur_pos[2]], VisitedType):
        return 0
    elif dbm[cur_pos[0]][cur_pos[1]][cur_pos[2]] != False:
        dbm[cur_pos[0]][cur_pos[1]][cur_pos[2]] = VisitedType(
            dbm[cur_pos[0]][cur_pos[1]][cur_pos[2]])
        return 1

    n_visits = 0
    # While the trail still has items in it...
    while len(trail) > 0:
        # If the current position has no more places to move to, pop it off (we're done with it),
        # move back to the previous position, and continue.
        if len(trail[-1][1]) == 0:
            trail.pop()
            if len(trail) > 0:
                cur_pos = trail[-1][0]
        # If there are still directions to try from the current position...
        else:
            # Pop a direction
            d = trail[-1][1].pop()
            next_pos = TUPLE_ADD(cur_pos, d)

            # Check to see if the next position is one of the following:
            # - Is out of bounds in any component (in which case we do nothing and iterate)
            if not TUPLE_LE(
                (-1, -1, -1), next_pos) or not TUPLE_LE(next_pos, M):
                pass
            # - Is empty, but has been visited (in which case we don't move there): None
            elif dbm[next_pos[0]][next_pos[1]][next_pos[2]] is None:
                pass
            # - Is occupied and has been visited (in which case we do nothing and iterate): (True,)
            elif isinstance(dbm[next_pos[0]][next_pos[1]][next_pos[2]],
                            VisitedType):
                pass
            # - Is occupied, but has not been visited (in which case we dont move there, and mark
            #   if as occupied and visited): True
            elif dbm[next_pos[0]][next_pos[1]][next_pos[2]] == False:
                # If the tests pass, this becomes the new current position, so push it, and a
                # fresh set of directions to try at this new position.
                trail.append((next_pos, copy(UNIT_VECTORS)))
                cur_pos = next_pos

                # If we reach this point, then the current position should be marked as having been
                # visited by the flood fill. This is distinguished from an unvisited location by
                # being set to None as opposed to False.
                #
                # This line will only execute when visiting valid new empty spaces.
                dbm[cur_pos[0]][cur_pos[1]][cur_pos[2]] = None
            else:
                n_visits += 1
                dbm[next_pos[0]][next_pos[1]][next_pos[2]] = VisitedType(
                    dbm[next_pos[0]][next_pos[1]][next_pos[2]])

    return n_visits


def flood_hollow_dbm(dbm, positions, test=lambda _: False):
    """
    Given a dense boolean matrix of which positions are filled, and the positions mapping, perform
    a flood fill from every point on the exterior of the bounding box, and remove all positions
    that aren't touched by the flood.
    """
    class __VisitedPosition:
        def __init__(self, val):
            self.val = val

    # For each positions on the edge of the bounding box, start to run a flood.
    # If the starting position is occupied, do nothing and move to the next.
    # If the starting position has been filled by a previous flood, do nothing and move on.
    x = len(dbm)
    y = len(dbm[0])
    z = len(dbm[0][0])

    starting_positions = [
        t for u in [[(0, j, k), (x - 1, j, k)] for j in range(y)
                    for k in range(z)] for t in u
    ]
    starting_positions += [
        t for u in [[(i, 0, k), (i, y - 1, k)] for i in range(x)
                    for k in range(z)] for t in u
    ]
    starting_positions += [
        t for u in [[(i, j, 0), (i, j, z - 1)] for i in range(x)
                    for j in range(y)] for t in u
    ]
    starting_positions = list(set(starting_positions))

    for start in starting_positions:
        iterative_flood_fill(dbm, start, __VisitedPosition)

    # Now that we've flooded the exterior, go through and prune anything that wasn't touched in the
    # flood. Also reset any point that was flooded but not originally filled to empty.
    pruned_positions = set()
    for i in range(x):
        for j in range(y):
            for k in range(z):
                # Visited locations that were occupied will be set to (True, <original_value>), and
                # this is exactly the set of locations that should be set to True; all others
                # should be set to False
                if isinstance(dbm[i][j][k], __VisitedPosition):
                    dbm[i][j][k] = dbm[i][j][k].val
                # The exception tot he above rule is the test parameter. If this evaluates to True
                # on the value of a given position, then that position should be spared.
                elif test(dbm[i][j][k]):
                    pass
                else:
                    # If the position is True, then it was set, but not visited and should be
                    # deleted.
                    if dbm[i][j][k] not in (None, False):
                        pruned_positions.add((k, j, i))
                    dbm[i][j][k] = False

    if isinstance(positions, list):
        return dbm, [
            pos for pos in positions if tuple(pos) not in pruned_positions
        ]
    elif isinstance(positions, dict):
        return dbm, dict([
            (pos, positions[pos]) for pos in
            [p for p in positions.keys() if p not in pruned_positions]
        ])


def csv_to_array(csv):
    return [[int(float(i)) for i in l.strip().split(",")]
            for l in csv.strip().split("\n")]
