# Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
# See the LICENSE.txt file at the top-level directory of this distribution.

import glob

env = Environment(LINK="gfortran", F90FLAGS="-O3 -ffast-math -g -Wall")
sources = glob.glob("*.f90")

env.Program("vd_test.x", sources)
