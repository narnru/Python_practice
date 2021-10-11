# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 20:57:57 2021

@author: Nikit
"""

import pygimli as pg
import pygimli.meshtools as mt
import time

world = mt.createWorld(start=[-20, 0], end=[20, -16], layers=[-2, -8],
                       worldMarker=False)
# Create a heterogeneous block
block = mt.createRectangle(start=[-6, -3.5], end=[6, -6.0],
                           marker=4,  boundaryMarker=10, area=0.1)
# Merge geometrical entities
geom = world + block
# pg.show(geom, boundaryMarker=True)

start = time.time()

mesh = mt.createMesh(geom, quality=33, area=0.01, smooth=[1, 10])
pg.show(mesh)

print(time.time() - start)

T = pg.solver.solveFiniteElements(mesh,
                                  a={1: 1.0, 2: 2.0, 3: 3.0, 4:0.1},
                                  bc={'Dirichlet': {8: 1.0, 4: 0.0}}, verbose=True)

ax, _ = pg.show(mesh, data=T, label='Temperature $T$', cMap="hot_r")
pg.show(geom, ax=ax, fillRegion=False)