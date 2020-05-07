#!/usr/bin/env python
"""
This script was generated as part of a driver_script file by the
setup_testcases.py script.
"""
import sys
import os
import shutil
import glob
import subprocess
import argparse


# This script was generated by setup_testcases.py as part of a driver_script
# file.
os.environ['PYTHONUNBUFFERED'] = '1'
parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--no_init_step1", dest="no_init_step1",
                    help="If set, init_step1 case will not be run during "
                         "execution of this script.",
                    action="store_true")
parser.add_argument("--finalize_init_step1", dest="finalize_init_step1",
                    help="If set, init_step1 case will have symlinks replaced "
                         "with the files they point to, this occurs after any "
                         "case runs that have been requested.",
                    action="store_true")
parser.add_argument("--no_forward1", dest="no_forward1",
                    help="If set, forward1 case will not be run during "
                         "execution of this script.",
                    action="store_true")
parser.add_argument("--finalize_forward1", dest="finalize_forward1",
                    help="If set, forward1 case will have symlinks replaced "
                         "with the files they point to, this occurs after any "
                         "case runs that have been requested.",
                    action="store_true")
parser.add_argument("--no_init_step2", dest="no_init_step2",
                    help="If set, init_step2 case will not be run during "
                         "execution of this script.",
                    action="store_true")
parser.add_argument("--finalize_init_step2", dest="finalize_init_step2",
                    help="If set, init_step2 case will have symlinks replaced "
                         "with the files they point to, this occurs after any "
                         "case runs that have been requested.",
                    action="store_true")
parser.add_argument("--no_forward2", dest="no_forward2",
                    help="If set, forward2 case will not be run during "
                         "execution of this script.",
                    action="store_true")
parser.add_argument("--finalize_forward2", dest="finalize_forward2",
                    help="If set, forward2 case will have symlinks replaced "
                         "with the files they point to, this occurs after any "
                         "case runs that have been requested.",
                    action="store_true")
parser.add_argument("--no_init_step3", dest="no_init_step3",
                    help="If set, init_step3 case will not be run during "
                         "execution of this script.",
                    action="store_true")
parser.add_argument("--finalize_init_step3", dest="finalize_init_step3",
                    help="If set, init_step3 case will have symlinks replaced "
                         "with the files they point to, this occurs after any "
                         "case runs that have been requested.",
                    action="store_true")
parser.add_argument("--no_forward3", dest="no_forward3",
                    help="If set, forward3 case will not be run during "
                         "execution of this script.",
                    action="store_true")
parser.add_argument("--finalize_forward3", dest="finalize_forward3",
                    help="If set, forward3 case will have symlinks replaced "
                         "with the files they point to, this occurs after any "
                         "case runs that have been requested.",
                    action="store_true")
parser.add_argument("--no_init_step4", dest="no_init_step4",
                    help="If set, init_step4 case will not be run during "
                         "execution of this script.",
                    action="store_true")
parser.add_argument("--finalize_init_step4", dest="finalize_init_step4",
                    help="If set, init_step4 case will have symlinks replaced "
                         "with the files they point to, this occurs after any "
                         "case runs that have been requested.",
                    action="store_true")
parser.add_argument("--no_forward4", dest="no_forward4",
                    help="If set, forward4 case will not be run during "
                         "execution of this script.",
                    action="store_true")
parser.add_argument("--finalize_forward4", dest="finalize_forward4",
                    help="If set, forward4 case will have symlinks replaced "
                         "with the files they point to, this occurs after any "
                         "case runs that have been requested.",
                    action="store_true")
parser.add_argument("--no_analysis", dest="no_analysis",
                    help="If set, analysis case will not be run during "
                         "execution of this script.",
                    action="store_true")
parser.add_argument("--finalize_analysis", dest="finalize_analysis",
                    help="If set, analysis case will have symlinks replaced "
                         "with the files they point to, this occurs after any "
                         "case runs that have been requested.",
                    action="store_true")

args = parser.parse_args()
base_path = os.getcwd()
dev_null = open('/dev/null', 'w')
error = False

if not args.no_init_step1:
    os.chdir(base_path)
    os.chdir('init_step1')

    print(" * Running init_step1")
    # Run command is:
    # ./run.py
    subprocess.check_call(['./run.py'], stdout=dev_null, stderr=None)
    print("     Complete")
if not args.no_forward1:
    os.chdir(base_path)
    os.chdir('forward1')

    print(" * Running forward1")
    # Run command is:
    # ./run.py
    subprocess.check_call(['./run.py'], stdout=dev_null, stderr=None)
    print("     Complete")
if not args.no_init_step2:
    os.chdir(base_path)
    os.chdir('init_step2')

    print(" * Running init_step2")
    # Run command is:
    # ./run.py
    subprocess.check_call(['./run.py'], stdout=dev_null, stderr=None)
    print("     Complete")
if not args.no_forward2:
    os.chdir(base_path)
    os.chdir('forward2')

    print(" * Running forward2")
    # Run command is:
    # ./run.py
    subprocess.check_call(['./run.py'], stdout=dev_null, stderr=None)
    print("     Complete")
if not args.no_init_step3:
    os.chdir(base_path)
    os.chdir('init_step3')

    print(" * Running init_step3")
    # Run command is:
    # ./run.py
    subprocess.check_call(['./run.py'], stdout=dev_null, stderr=None)
    print("     Complete")
if not args.no_forward3:
    os.chdir(base_path)
    os.chdir('forward3')

    print(" * Running forward3")
    # Run command is:
    # ./run.py
    subprocess.check_call(['./run.py'], stdout=dev_null, stderr=None)
    print("     Complete")
if not args.no_init_step4:
    os.chdir(base_path)
    os.chdir('init_step4')

    print(" * Running init_step4")
    # Run command is:
    # ./run.py
    subprocess.check_call(['./run.py'], stdout=dev_null, stderr=None)
    print("     Complete")
if not args.no_forward4:
    os.chdir(base_path)
    os.chdir('forward4')

    print(" * Running forward4")
    # Run command is:
    # ./run.py
    subprocess.check_call(['./run.py'], stdout=dev_null, stderr=None)
    print("     Complete")
if not args.no_analysis:
    os.chdir(base_path)
    os.chdir('analysis')

    print(" * Running analysis")
    # Run command is:
    # ./run.py
    subprocess.check_call(['./run.py'], stdout=dev_null, stderr=None)
    print("     Complete")
if error:
    sys.exit(1)
if args.finalize_init_step1:
    old_dir = os.getcwd()
    os.chdir("init_step1")
    file_list = glob.glob("*")
    for file in file_list:
        if os.path.islink(file):
            link_path = os.readlink(file)
            os.unlink(file)
            shutil.copyfile(link_path, file)
    os.chdir(old_dir)

if args.finalize_forward1:
    old_dir = os.getcwd()
    os.chdir("forward1")
    file_list = glob.glob("*")
    for file in file_list:
        if os.path.islink(file):
            link_path = os.readlink(file)
            os.unlink(file)
            shutil.copyfile(link_path, file)
    os.chdir(old_dir)

if args.finalize_init_step2:
    old_dir = os.getcwd()
    os.chdir("init_step2")
    file_list = glob.glob("*")
    for file in file_list:
        if os.path.islink(file):
            link_path = os.readlink(file)
            os.unlink(file)
            shutil.copyfile(link_path, file)
    os.chdir(old_dir)

if args.finalize_forward2:
    old_dir = os.getcwd()
    os.chdir("forward2")
    file_list = glob.glob("*")
    for file in file_list:
        if os.path.islink(file):
            link_path = os.readlink(file)
            os.unlink(file)
            shutil.copyfile(link_path, file)
    os.chdir(old_dir)

if args.finalize_init_step3:
    old_dir = os.getcwd()
    os.chdir("init_step3")
    file_list = glob.glob("*")
    for file in file_list:
        if os.path.islink(file):
            link_path = os.readlink(file)
            os.unlink(file)
            shutil.copyfile(link_path, file)
    os.chdir(old_dir)

if args.finalize_forward3:
    old_dir = os.getcwd()
    os.chdir("forward3")
    file_list = glob.glob("*")
    for file in file_list:
        if os.path.islink(file):
            link_path = os.readlink(file)
            os.unlink(file)
            shutil.copyfile(link_path, file)
    os.chdir(old_dir)

if args.finalize_init_step4:
    old_dir = os.getcwd()
    os.chdir("init_step4")
    file_list = glob.glob("*")
    for file in file_list:
        if os.path.islink(file):
            link_path = os.readlink(file)
            os.unlink(file)
            shutil.copyfile(link_path, file)
    os.chdir(old_dir)

if args.finalize_forward4:
    old_dir = os.getcwd()
    os.chdir("forward4")
    file_list = glob.glob("*")
    for file in file_list:
        if os.path.islink(file):
            link_path = os.readlink(file)
            os.unlink(file)
            shutil.copyfile(link_path, file)
    os.chdir(old_dir)

if args.finalize_analysis:
    old_dir = os.getcwd()
    os.chdir("analysis")
    file_list = glob.glob("*")
    for file in file_list:
        if os.path.islink(file):
            link_path = os.readlink(file)
            os.unlink(file)
            shutil.copyfile(link_path, file)
    os.chdir(old_dir)

sys.exit(0)
