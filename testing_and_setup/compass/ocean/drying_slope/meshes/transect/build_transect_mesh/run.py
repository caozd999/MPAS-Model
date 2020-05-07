#!/usr/bin/env python

# This script was generated from setup_testcases.py as part of a config file

import sys
import os
import shutil
import glob
import subprocess


dev_null = open('/dev/null', 'w')

# Run command is:
# planar_hex --nx 4 --ny 116 --dc 80.0 --npy -o grid.nc
subprocess.check_call(['planar_hex', '--nx', '4', '--ny', '116', '--dc',
                       '80.0', '--npy', '-o', 'grid.nc'])

# Run command is:
# MpasCellCuller.x grid.nc culled_mesh.nc
subprocess.check_call(['MpasCellCuller.x', 'grid.nc', 'culled_mesh.nc'])

# Run command is:
# MpasMeshConverter.x culled_mesh.nc transectmesh.nc
subprocess.check_call(['MpasMeshConverter.x', 'culled_mesh.nc',
                       'transectmesh.nc'])
