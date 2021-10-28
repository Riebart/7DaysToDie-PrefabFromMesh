#!/usr/bin/env python3

import argparse
from multiprocessing import set_executable
import struct
import random
import sys
import time
import math
from stl_to_pts import *

# from _typeshed import FileDescriptor

# Predetermined lambdas to remap dimensions, based on the 6 possible options:
DIM_REMAP_LAMBDA = {
    "123": lambda p: p,
    "132": lambda p: (p[0], p[2], p[1]),
    "213": lambda p: (p[1], p[0], p[2]),
    "231": lambda p: (p[1], p[2], p[0]),
    "312": lambda p: (p[2], p[0], p[1]),
    "321": lambda p: (p[2], p[1], p[0])
}


class TTSFile(object):
    """
    Parse a TTS file given a file descriptor, allowing operations to take place on high level constructs, and write the resulting output back to a byte array, or file descriptor.
    """
    @staticmethod
    def block_id_flags(block_value):
        return (block_value & 0x00007fff, block_value & 0xffff8000)

    @staticmethod
    def block_value(block_id, block_flags):
        return block_id + (block_flags << 15)

    def __init__(self, tts_fd=None):
        magic_header = tts_fd.read(4).decode("utf-8")
        assert magic_header == "tts\0"
        (self.version, self.x_size, self.y_size,
         self.z_size) = struct.unpack("<LHHH", fd.read(10))

        total_volume_size = self.x_size * self.y_size * self.z_size

        block_data = struct.unpack("<" + ("L" * total_volume_size),
                                   fd.read(total_volume_size * 4))

        self.blocks = [[[
            TTSFile.block_id_flags(block_data[z * self.x_size * self.y_size +
                                              y * self.x_size + x])
            for x in range(self.x_size)
        ] for y in range(self.y_size)] for z in range(self.z_size)]

        # The way we're assuming the ending metadata works is overly simplistic
        self.per_block_end_headers = list()
        self.per_block_end_headers.append(fd.read(1))
        fd.read(total_volume_size - 1)
        self.per_block_end_headers.append(fd.read(1))
        fd.read(total_volume_size - 1)
        self.per_block_end_headers.append(fd.read(1))
        fd.read(total_volume_size - 1)

        # The final footer is just a pile of random goop that seems to be:
        # - A 32-bit integer that is the total volume in blocks divided by 8, rounded up
        # - N bytes, where N is the previous total volume divided by 8 rounded up
        # - 2 more zero bytes at the end
        #self.footer = fd.read()

        self.observed_blocks = list(
            set([
                self.blocks[z][y][x] for z in range(self.z_size)
                for y in range(self.y_size) for x in range(self.x_size)
                if self.blocks[z][y][x] != (0, 0)
            ]))

    def write(self, output_fd):
        fd.write(b"tts\0")
        fd.write(struct.pack("<L", self.version))
        fd.write(struct.pack("<H", self.x_size))
        fd.write(struct.pack("<H", self.y_size))
        fd.write(struct.pack("<H", self.z_size))

        for z in self.blocks:
            for y in z:
                for x in y:
                    fd.write(struct.pack("<L", TTSFile.block_value(*x)))

        total_blocks = self.x_size * self.y_size * self.z_size

        for meta_val in self.per_block_end_headers:
            fd.write(meta_val * total_blocks)

        fd.write(struct.pack("<L", math.ceil(total_blocks / 8.0)))
        fd.write(b"\0" * math.ceil(total_blocks / 8.0))

        fd.write(b"\0\0")
        fd.flush()

    def clear(self):
        """
        Reset all block values to air
        """
        self.blocks = [[[(0, 0) for _ in range(self.x_size)]
                        for _ in range(self.y_size)]
                       for _ in range(self.z_size)]
        # for z in self.blocks:
        #     for y in z:
        #         for x in y:
        #             x = (0, 0)

    def resize(self,
               x_size: int,
               y_size: int,
               z_size: int,
               x_offset: int = 0,
               y_offset: int = 0,
               z_offset: int = 0):
        """
        Resize a TTS dense block cloud
        """
        new_blocks = [[[(0, 0) for _ in range(x_size)] for _ in range(y_size)]
                      for _ in range(z_size)]

        # Copy the original data into the new block cloud
        for zi in range(min(self.z_size, z_size - z_offset)):
            for yi in range(min(self.y_size, y_size - y_offset)):
                for xi in range(min(self.x_size, x_size - x_offset)):
                    new_blocks[zi + z_offset][yi + y_offset][
                        xi + x_offset] = self.blocks[zi][yi][xi]

        self.blocks = new_blocks
        self.x_size = x_size
        self.y_size = y_size
        self.z_size = z_size

    def paint(self, positions):
        for p in positions:
            self.blocks[p[2]][p[1]][p[0]] = self.observed_blocks[0]

    def __str__(self):
        return str((self.version, self.x_size, self.y_size, self.z_size,
                    len(self.blocks), self.blocks[0][0][0],
                    self.per_block_end_headers, self.footer))

    def draw(self, fname="output.png"):
        from os import environ
        environ['PYGAME_HIDE_SUPPORT_PROMPT'] = '1'
        import pygame

        colors = {}
        block_size = (10, 10, 3)

        image_size_x = ((1 * self.x_size * block_size[0]) +
                        (self.z_size * block_size[2]) + block_size[0])
        image_size_y = ((1 * self.y_size * block_size[1]) +
                        (self.z_size * block_size[2]) + block_size[1])
        image = pygame.surface.Surface((image_size_x, image_size_y))

        z = 0
        for each_layer in self.blocks:
            y = 0
            for each_row in each_layer:
                x = 0
                for each_block in each_row:
                    if each_block != 0 and (isinstance(each_block, tuple)
                                            and each_block[0] != 0):
                        if each_block not in colors:
                            colors[each_block] = (random.randint(0, 255),
                                                  random.randint(0, 255),
                                                  random.randint(0, 255))
                        draw_color = colors[each_block]
                        draw_x = (x * block_size[0]) + (z * block_size[2])
                        draw_y = (-y * block_size[1]) + (
                            z * block_size[2]) + self.y_size * block_size[1]
                        draw_rect = (draw_x, draw_y, block_size[0],
                                     block_size[1])

                        pygame.draw.rect(image, draw_color, draw_rect, 0)
                        pygame.draw.rect(image, (0, 0, 0), draw_rect, 1)
                    x += 1
                y += 1
            z += 1

        print("Block ids and colors: " + str(colors))
        pygame.image.save(image, fname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stl-file", help="Path to STL input model")
    parser.add_argument("--tts-out-file", help="Path to output TTS file")
    parser.add_argument(
        "--prefab-size",
        type=int,
        help=
        "Size, in voxels, of the resulting prefab, applied to the longest dimension of the input model file",
        required=True)
    parser.add_argument(
        "--dimension-remap",
        help=
        "Some permutation of 123, indicating that dimensions of points should be remapped to rotate the object appropriately.",
        required=False,
        default="123")
    pargs = parser.parse_args()

    with open(pargs.tts_out_file, "rb") as fd:
        tts = TTSFile(fd)
        tts.draw("outputA.png")

    with open(pargs.stl_file, "rb") as fd:
        t0 = time.time()
        triangles = STLFile.read_triangles(fd)

    sys.stderr.write("Reading model took %s seconds.\n" %
                     str(time.time() - t0))
    sys.stderr.write("Model has %d triangles\n" % len(triangles))

    bounds = triangle_list_bounds(triangles)
    sys.stderr.write("Model bounds: %s\n" % str(bounds))

    # To assist with ensuring symmetry, shift the points so that the centroid
    # of the model is at the origin. Find the midpoint along each of dimensions
    # of the cube spanned by the bounds of the model, and subtract that midpoint
    # from each triangle coordinate.
    origin_offset = [-sum(b) / 2 for b in bounds]
    for t in triangles:
        t.shift(origin_offset)

    # For clarity, show the transposed model bounds, which should be symmetric.
    bounds = triangle_list_bounds(triangles)
    sys.stderr.write("Translated model bounds: %s\n" % str(bounds))

    longest_dim = max([i[1] - i[0] for i in bounds])
    resolution = longest_dim / (abs(pargs.prefab_size) - 1)
    sys.stderr.write("Computed spatial resolution in model-space: %f\n" %
                     resolution)

    sys.stderr.write("Splitting triangles...\n")
    t0 = time.time()
    pts = split_tris(triangles, resolution)
    sys.stderr.write("Triangle to point refinement took %s seconds.\n" %
                     str(time.time() - t0))
    sys.stderr.write("Split %d triangles into %d points.\n" %
                     (len(triangles), len(pts)))

    # Remap the dimensions of the STL file to orient it as specified.
    pts = list(map(DIM_REMAP_LAMBDA[pargs.dimension_remap], pts))

    # Step 1: Figure out how big the bounding box is, calculate the two opposing corners.
    m, M = bounding_box(pts)

    # Step 2: Move all positions by the minimal corner, putting one corner at the origin,
    # and move the maximal corner
    length, width, height = list_subtract(M, m)
    length += 1
    width += 1
    height += 1
    sys.stderr.write("Dimensions of resulting prefab: %d %d %d\n" %
                     (length, width, height))
    sys.stderr.write("Resulting prefab chunk count: %d" %
                     (math.ceil(length * width * height / 8.0)))
    ptsT = [list_subtract(p, m) for p in pts]

    tts.resize(length, width, height)
    tts.draw("outputB.png")

    tts.clear()
    tts.paint(ptsT)
    tts.draw("outputC.png")

    with open("axis-test.tts", "wb") as fd:
        tts.write(fd)
