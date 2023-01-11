#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a program to convert the graph data from
Stanford’s Large Network Dataset Collection to LARC
Matrix Market mmio format.
"""

import argparse
import configparser
import sys
import os
import os.path
import re

__author__ = "IDA/CCS LARC Team"
__copyright__ = "Copyright (c) 2021 by Institute for Defense Analyses.  All rights reserved."
__version__ = "1.0.0"
__date__ = "March 9, 2021"


# Global verbose flag.
# If something other than None is specified, it will override any command line arguments.
verbose = None


def do_one_file(infp, outfp, infilename, outfilename, profile):
    print("%%MatrixMarket matrix coordinate integer {0}".format(profile["mmio_type"]), file=outfp)
    print("% Generated from input file: {0}".format(infilename), file=outfp)

    num_edges = 0
    max_vertex = 0

    pattern_two_nums = re.compile("([0-9]+)(?:[ \t]*,[ \t]*|[ \t]+)([0-9]+)")

    for linenumber, line in enumerate(infp):
        stripline = line.strip()
        if stripline == "":
            pass
        elif stripline[0] == '#':
            print("% [{0}] >>> {1}".format(linenumber+1, stripline), file=outfp)
        elif stripline[0] == '%':
            print("% [{0}] >>> {1}".format(linenumber+1, stripline), file=outfp)
        else:
            match_result = pattern_two_nums.fullmatch(stripline)
            if match_result:
                first_entry = int(match_result.group(1))
                second_entry = int(match_result.group(2))
                first_entry += profile["vertex_shift"]
                second_entry += profile["vertex_shift"]
                max_vertex = max(max_vertex, first_entry, second_entry)

                if first_entry < second_entry:
                    edge_proc = profile["up_edges"]
                elif first_entry > second_entry:
                    edge_proc = profile["down_edges"]
                else:
                    edge_proc = profile["self_edges"]

                if edge_proc in ("Output", "Reverse"):
                    num_edges += 1
                elif edge_proc in ("Double"):
                    num_edges += 2
                elif edge_proc == "Skip":
                    print("% [{0}] SKIP >>> {1}".format(linenumber+1, stripline), file=outfp)
                elif edge_proc == "Fail":
                    assert False, "Bad edge {0} {1} encountered; translation failed.".format(
                            first_entry, second_entry )
            else:
                print("% [{0}] ??? >>> {1}".format(linenumber+1, stripline), file=outfp)

    print("{0} {0} {1}".format(max_vertex, num_edges), file=outfp)

    infp.seek(0)
    for line in infp:
        stripline = line.strip()
        if stripline == "":
            pass
        elif stripline[0] == '#':
            pass
        elif stripline[0] == '%':
            pass
        else:
            match_result = pattern_two_nums.fullmatch(stripline)
            if match_result:
                first_entry = int(match_result.group(1))
                second_entry = int(match_result.group(2))
                first_entry += profile["vertex_shift"]
                second_entry += profile["vertex_shift"]

                if first_entry < second_entry:
                    edge_proc = profile["up_edges"]
                elif first_entry > second_entry:
                    edge_proc = profile["down_edges"]
                else:
                    edge_proc = profile["self_edges"]

                if edge_proc == "Output":
                    print("{0} {1} 1".format(first_entry, second_entry), file=outfp)
                elif edge_proc == "Reverse":
                    print("{1} {0} 1".format(first_entry, second_entry), file=outfp)
                elif edge_proc in ("Double"):
                    print("{0} {1} 1".format(first_entry, second_entry), file=outfp)
                    print("{1} {0} 1".format(first_entry, second_entry), file=outfp)
                elif edge_proc == "Skip":
                    pass
                elif edge_proc == "Fail":
                    assert False
            else:
                pass


def scan_input_file(infp, infilename):
    print("SCANNING Input file: {0}".format(infilename))

    num_blank_lines = 0
    num_comment_lines_sharp = 0
    num_comment_lines_percent = 0
    num_bad_lines = 0
    num_edge_lines = 0
    num_edge_lines_increasing = 0
    num_edge_lines_decreasing = 0
    num_edge_lines_equal = 0
    vertex_set = set()

    pattern_two_nums = re.compile("([0-9]+)(?:[ \t]*,[ \t]*|[ \t]+)([0-9]+)")

    for linenumber, line in enumerate(infp):
        stripline = line.strip()
        if stripline == "":
            num_blank_lines += 1
        elif stripline[0] == '#':
            print("[{0}] >>> {1}".format(linenumber+1, stripline))
            num_comment_lines_sharp += 1
        elif stripline[0] == '%':
            print("[{0}] >>> {1}".format(linenumber+1, stripline))
            num_comment_lines_percent += 1
        else:
            match_result = pattern_two_nums.fullmatch(stripline)
            if match_result:
                first_entry = int(match_result.group(1))
                second_entry = int(match_result.group(2))
                vertex_set.add(first_entry)
                vertex_set.add(second_entry)
                if first_entry == second_entry:
                    print("[{0}] SELF >>> {1}".format(linenumber+1, stripline))
                    num_edge_lines_equal += 1
                elif first_entry < second_entry:
                    if num_edge_lines_increasing < 3:
                        print("[{0}] INCREASING >>> {1}".format(linenumber+1, stripline))
                    num_edge_lines_increasing += 1
                else:
                    if num_edge_lines_decreasing < 3:
                        print("[{0}] DECREASING >>> {1}".format(linenumber+1, stripline))
                    num_edge_lines_decreasing += 1
                num_edge_lines += 1
            else:
                print("[{0}] BAD >>> {1}".format(linenumber+1, stripline))
                num_bad_lines += 1

    print()
    print("# of blank lines = {0}.".format(num_blank_lines))
    print("# of (sharp) comment lines = {0}.".format(num_comment_lines_sharp))
    print("# of (percent) comment lines = {0}.".format(num_comment_lines_percent))
    print("# of bad lines = {0}.".format(num_bad_lines))
    print("# of edge lines = {0}, up = {1}, down = {2}, self loop = {3}.".format(
            num_edge_lines, num_edge_lines_increasing, num_edge_lines_decreasing, num_edge_lines_equal ))
    print("# of distinct vertices = {0}, min = {1}, max = {2}.".format(
            len(vertex_set), min(vertex_set), max(vertex_set) ))
    print("SCAN COMPLETE!")


def create_default_configuration_file(profile_filename):
    contents = \
"""[General]
version = 1.0.0

[Translate]
# Allowed values for mmio_type are 'general' and 'symmetric'.
# When 'symmetric' is used, only the data for the lower triangular
# portion should be generated, i.e. for row index >= column index,
# where the row index is given first and the column index second.
mmio_type = general

# The vertex_shift will be added to the vertex numbers when they are read in.
# After the addition, the result should be an index starting at 1.
vertex_shift = 0

# Supported actions are: 'Output', 'Reverse', 'Double', 'Skip', 'Fail'.
up_edges = Output
down_edges = Output
self_edges = Output
"""
    try:
        with open(profile_filename, "w") as cfp:
           cfp.write(contents)
    except Exception as err:
        print("ERROR: Exception encountered while trying to create default configuration file '{0}':".format(
                profile_filename ))
        print(" ---> {0}".format(err))
        print("PROGRAM TERMINATING!")
        sys.exit(1)
    

def process_configuration_profile(profile_filename):
    config = configparser.ConfigParser()
    config.read(profile_filename)
    if "Translate" in config:
        trans = config["Translate"]
    else:
        print("ERROR: Configuration file is missing a [Translate] section.")
        print("PROGRAM TERMINATING!")
        sys.exit(1)

    profile = {}
    profile["mmio_type"] = trans.get("mmio_type")
    profile["vertex_shift"] = trans.getint("vertex_shift")
    profile["up_edges"] = trans.get("up_edges")
    profile["down_edges"] = trans.get("down_edges")
    profile["self_edges"] = trans.get("self_edges")

    if profile["mmio_type"] == None:
        print("ERROR: Configuration file does not contain entry for 'mmio_type' under [Translate].")
        print("PROGRAM TERMINATING!")
        sys.exit(1)
    elif profile["mmio_type"] not in ("general", "symmetric"):
        print("ERROR: Invalid entry '{0}' for 'mmio_type' in configuration file.".format(
                profile["mmio_type"] ))
        print("PROGRAM TERMINATING!")
        sys.exit(1)

    if profile["vertex_shift"] == None:
        profile["vertex_shift"] = 0

    edge_processing_options = ("Output", "Reverse", "Double", "Skip", "Fail")

    if profile["up_edges"] == None:
        print("ERROR: Configuration file does not contain entry for 'up_edges' under [Translate].")
        print("PROGRAM TERMINATING!")
        sys.exit(1)
    elif profile["up_edges"] not in edge_processing_options:
        print("ERROR: Invalid entry '{0}' for 'up_edges' in configuration file.".format(
                profile["up_edges"] ))
        print("PROGRAM TERMINATING!")
        sys.exit(1)

    if profile["down_edges"] == None:
        print("ERROR: Configuration file does not contain entry for 'down_edges' under [Translate].")
        print("PROGRAM TERMINATING!")
        sys.exit(1)
    elif profile["down_edges"] not in edge_processing_options:
        print("ERROR: Invalid entry '{0}' for 'down_edges' in configuration file.".format(
                profile["down_edges"] ))
        print("PROGRAM TERMINATING!")
        sys.exit(1)

    if profile["self_edges"] == None:
        print("ERROR: Configuration file does not contain entry for 'self_edges' under [Translate].")
        print("PROGRAM TERMINATING!")
        sys.exit(1)
    elif profile["self_edges"] not in edge_processing_options:
        print("ERROR: Invalid entry '{0}' for 'self_edges' in configuration file.".format(
                profile["self_edges"] ))
        print("PROGRAM TERMINATING!")
        sys.exit(1)

    return profile


def execute_command(argstruct):
    infilename = argstruct.infile
    outfilename = argstruct.outfile
    try:
        infp = open(infilename, "r")
    except Exception as err:
        print("ERROR: Exception encountered while trying to open input file '{0}':".format(infilename))
        print(" ---> {0}".format(err))
        print("PROGRAM TERMINATING!")
    else:
        with infp:
            if argstruct.scan:
                scan_input_file(infp, infilename)
            else:
                try:
                    outfp = open(outfilename, "w")
                except Exception as err:
                    print("ERROR: Exception encountered while trying to open output file '{0}':".format(
                            outfilename ))
                    print(" ---> {0}".format(err))
                    print("PROGRAM TERMINATING!")
                else:
                    with outfp:
                        outfiledir = os.path.dirname(outfilename)
                        profile_filename = os.path.join(outfiledir, "convert_profile.ini")
                        if not os.path.exists(profile_filename):
                            if verbose > 0:
                                print("WARNING: Creating default configuration file '{0}'.".format(profile_filename))
                            create_default_configuration_file(profile_filename)
                        profile = process_configuration_profile(profile_filename)
                        do_one_file(infp, outfp, infilename, outfilename, profile)


def expand_and_validate_args(argstruct):
    # Set global verbose flag.
    global verbose
    if verbose == None:
        verbose = argstruct.verbose
        if argstruct.quiet:
            if argstruct.verbose == 1:
                verbose = 0
            else:
                print("WARNING: Ignoring '--quiet' flag due to '--verbose' specification on command line.")
    elif verbose > 0:
        print("INFO: Ignoring command line quiet and/or verbose settings.")
        print(" ---> Using hard-coded verbose setting {0} instead.".format(verbose))

    # Nothing else to do if user specified a scan only.
    if argstruct.scan:
        return

    # Supply default outfilename if user did not supply one.
    if argstruct.outfile == None:
        argstruct.defaulted_outfile = True
        suffix_position = argstruct.infile.rfind('.')
        if suffix_position < 0:
            argstruct.outfile = argstruct.infile + ".mmio"
        else:
            argstruct.outfile = argstruct.infile[:suffix_position] + ".mmio"
    else:
        argstruct.defaulted_outfile = False

    # Check if outfile already exists.
    if os.path.exists(argstruct.outfile):
        if argstruct.force:
            if verbose > 0:
                print("WARNING: Overwriting output file '{0}'.".format(argstruct.outfile))
        else:
            print("ERROR: Output file '{0}' already exists.".format(argstruct.outfile))
            print("Specify '--force' to overwrite output file.")
            print("PROGRAM TERMINATING!")
            sys.exit(1)


def run_from_main():
    parser = argparse.ArgumentParser(description="Convert Stanford SNAP data file to LARC mmio.")
    parser.add_argument("infile", help="input file name")
    parser.add_argument("outfile", nargs="?", help="output file name")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing output file")
    parser.add_argument("-q", "--quiet", action="store_true", help="supress informational messages")
    parser.add_argument("-v", "--verbose", action="count", default=1, help="generate more informational messages")
    parser.add_argument("-s", "--scan", action="store_true", help="perform input file scan only")
    parser.add_argument("--version", action="version", version=__version__)
    argstruct = parser.parse_args()
    expand_and_validate_args(argstruct)
    execute_command(argstruct)


if __name__ == "__main__":
    run_from_main()

