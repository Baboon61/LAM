# A script that manage to create an horizontal bar visualisation of double
# strand break in LAM-HTGTS method

##############################IMPORTS##############################

import os
import sys
import getopt
import os.path
import re
import csv
import pandas as pd
import numpy as np
import math
import operator
import collections
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.ticker import MultipleLocator

# USED
#(0.0, 0.0, 0.0, 1.0) = black = intra
#(1.0, 0.0, 0.0, 1.0)) = red = inter


def usage():
    print('Usage:\n')
    print('\tpython ' +
          sys.argv[0] + ' -m <metadata file> -g <genome type> -p <postprocess directory> -o <output_mark> -e <method type> [-i <mark input>]')
    print('\t\t-h or --help : display this help')
    print('\t\t-m or --file_metadata : metadata file')
    print('\t\t-g or --genome : only filter librairies results with this genome')
    print('\t\t-p or --dir_post : postprocess directory')
    print('\t\t-o or --output_mark : mark added to the output file')
    print('\t\t-e (0,1) or --method (0,1) : 0 = output picture using counts, 1 = output picture using log2(counts)')
    print('\t\t-i or --input_mark : mark from input file')

##############################FUNCTIONS##############################

# ALLOW TO DIVIDE ALL VALUES WITH MAX_DISTANCE_COUNT_THEN REMOVE
# MAX_DISTANCE_COUNT AND MAX_DISTANCE FROM DICTIONNARY


def manageArray(dictionnary, max_distance_count, method):
    # print(dictionnary)
    del dictionnary['max_distance']
    if dictionnary['max_distance_count'] != 0:
        for key, value in dictionnary.items():
            if method == 0:
                # COUNTS
                dictionnary[key] = float(
                    float(dictionnary[key]) / float(max_distance_count))
            elif method == 1:
                # LOG2
                dictionnary[key] = float(float(
                    np.log2(dictionnary[key] + 1)) / float(np.log2(max_distance_count + 1)))
            else:
                print("Wrong method use 0 or 1")
                sys.exit(2)
    del dictionnary['max_distance_count']
    return dictionnary


def main(argv):

    file_metadata = ""
    dir_post = ""
    output_mark = ""
    genome = ""
    method = -1
    file_input_extension = ""
    file_output_extension = ""
    title_label = ""

    good_directory_list = []
    check_result_tlx = False
    check_result_filter_stats_txt = False
    total_junction = 0

    #pd.options.mode.chained_assignment = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'm:g:p:o:e:i:', [
                                   'file_metadata=', 'genome=', 'dir_post=', 'output_mark=', 'method=', 'input_mark=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

##############################OPTIONS##############################

    print("\n")
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-m', '--file_metadata'):
            file_metadata = arg
        elif opt in ('-g', '--genome'):
            genome = arg
        elif opt in ('-p', '--dir_post'):
            dir_post = arg
        elif opt in ('-o', '--output_mark'):
            output_mark = arg
        elif opt in ('-e', '--method'):
            method = arg
        elif opt in ('-i', '--input_mark'):
            input_mark = arg
        else:
            print("Error : Bad option -> " + opt)
            usage()
            sys.exit(2)

##############################CHECK UP/SET UP##############################

    # CHECK METADATA FILE
    if file_metadata == "" or not os.path.exists(file_metadata):
        print("Error : You have to set a metadata file !\n")
        usage()
        sys.exit(2)
    else:
        # READ METADATA FILE
        metadata = pd.read_table(file_metadata, sep='\t')

    # FILTER METADATA FILE IF GENOME INPUT
    if genome != "":
        metadata = metadata.loc[metadata['Assembly'] == genome]
        if metadata.empty:
            print("Error : This assembly does not exist in the metadata file !\n")
            usage()
            sys.exit(2)
    else:
        print("Error : You have to set a genome name !\n")
        usage()
        sys.exit(2)

    # CHECK POSTPROCESS DIRECTORY
    if not os.path.exists(dir_post):
        print("Error : You have to set a postprocess directory !\n")
        usage()
        sys.exit(2)
    else:
        if dir_post[-1] != "/":
            dir_post += "/"

    # C HECK METHOD
    try:
        method = int(method)
        if method == 0:
            title_label = "Number of double stranded breaks after the primer at each position divided by the maximum of breaks by library"
        elif method == 1:
            title_label = "Log2 of number of double stranded breaks after the primer at each position divided by log2 of the maximum of breaks by library"
        else:
            print("Error : Method option needs to be 0 or 1 !\n")
            usage()
            sys.exit(2)
    except:
        print("Error : You have to set an integer (0 or 1) to method option !\n")
        usage()
        sys.exit(2)

    # CHECK INPUT MARKS HISTORY
    if input_mark == "":
        print("Warning : You will process the raw file !\n")

    # CHECK OUTPUT MARK
    if output_mark == "":
        print("Error : You have to set an output mark !\n")
        usage()
        sys.exit(2)

    # SELECT INPUT FILES
    if input_mark == "":
        file_input_extension = ".tlx"
    else:
        file_input_extension = "_" + "_".join(input_mark.split(",")) + ".tlx"

    # TEST IF INPUT EXIST IN AT LEAST ON LIBRARY
    check_input_mark = False
    for library in metadata['Library'].tolist():
        if os.path.exists(dir_post + library + "/" + library + "_Legitimate" + file_input_extension) or os.path.exists(dir_post + library + "/" + library + "_Illegitimate" + file_input_extension):
            check_input_mark = True
    if not check_input_mark:
        print("Error : Your input marks can not localize a good legitimate or illegitimate file !\n")
        usage()
        sys.exit(2)

    # SELECT OUTPUT FILES
    if file_input_extension != "":
        file_output_extension = file_input_extension[
            :-4] + "_" + output_mark + ".png"
    else:
        file_output_extension = "_" + output_mark + ".png"

##############################PRINTS##############################

    print('\n-----------------------------------------')
    print('Metadata file : ' + file_metadata)
    print('Genome : ' + genome)
    print('Postprocess directory : ' + dir_post)
    print('Method : ' + str(method))
    print('Input file extension: ' + file_input_extension)
    print('Output file extension : ' + file_output_extension)
    print('-----------------------------------------\n')

##############################PROGRAMS##############################

    fig = plt.figure()
    ax = fig.add_subplot(111)
    d = {}
    yticks = []
    yticklabels = []
    track_number = 0
    yrange = []
    coff_leg_ille_graph = 0

    # LOOP OVER EACH LIBRARIES
    for library in metadata['Library'].tolist():
        # print(library)
        # CHECK DIRECTORY EXISTS
        if not os.path.exists(dir_post + library):
            print("Warning : " + dir_post +
                  " does not contains {" + library + "}")
            print("Warning :  {" + library + "} will not be filtered")
        else:
            # CHECK INPUT FILE EXISTS
            if os.path.exists(dir_post + library + "/" + library + "_Legitimate" + file_input_extension):
                df_legitimate = pd.read_csv(dir_post + library + "/" + library +
                                            "_Legitimate" + file_input_extension, sep='\t', header=0, index_col=None)
                df_legitimate = df_legitimate.drop(
                    columns=df_legitimate.columns[11:])

            else:
                print("Error : The Legitimate file for " +
                      library + " is missing !\n")
                usage()
                sys.exit(2)

            if os.path.exists(dir_post + library + "/" + library + "_Illegitimate" + file_input_extension):
                df_illegitimate = pd.read_csv(
                    dir_post + library + "/" + library + "_Illegitimate" + file_input_extension, sep='\t', header=0, index_col=None)
                df_illegitimate = df_illegitimate.drop(
                    columns=df_illegitimate.columns[11:])

            else:
                print("Error : The Illegitimate file for " +
                      library + " is missing !\n")
                usage()
                sys.exit(2)

    metadata = metadata.sort_values(['Library'], ascending=[True])
    global_list = dir_post + metadata['Library']
    # for i in global_list:
    #	print(i)

    nb_files = len(global_list)
    # DICTIONNARY WITH ILLEGITIMATES AND LEGITIMATES JUNCTIONS
    distance_dict = {"illegitimates": {}, "legitimates": {}}

    # OPEN TLX FILES ONE BY ONE
    for result_file in global_list:
        label = result_file.split("/")[-1]

        # YRANGE USED AT THE END TO HORIZONTAL SPACE BAR

        yrange.append(((track_number * 1.8) + (track_number *
                                               1.8) + 1.8 + coff_leg_ille_graph, 1))
        yrange.append(((track_number * 1.8) + (track_number *
                                               1.8) + 3.6 + coff_leg_ille_graph, 1))

        coff_leg_ille_graph += 1

        # DICTIONNARY LABEL INSIDE ILLEGITIMATES AND LEGITIMATES DICTIONNARIES
        distance_dict["illegitimates"][label] = {}
        distance_dict["legitimates"][label] = {}

        # with open(result_file+"/"+label+"_result.tlx", 'r') as f_tlx:
        with open(result_file + "/" + label + "_Legitimate" + file_input_extension, 'r') as f_tlx:
            f_tlx.readline()
            for line in f_tlx:
                if pd.isnull(metadata.loc[metadata['Library'] == label]['MID'].values[0]):
                    len_mid = 0
                else:
                    len_mid = len(metadata.loc[metadata['Library'] == label]['MID'].values[0])
                if pd.isnull(metadata.loc[metadata['Library'] == label]['Primer'].values[0]):
                    len_primer = 0
                else:
                    len_primer = len(metadata.loc[metadata['Library'] == label]['Primer'].values[0])
                distance = int(
                    int(line.split("\t")[9]) - int(line.split("\t")[8])) - (len_mid + len_primer)
                if distance not in distance_dict["legitimates"][label]:
                    distance_dict["legitimates"][label][distance] = 1
                else:
                    distance_dict["legitimates"][label][distance] += 1
                total_junction += 1

        with open(result_file + "/" + label + "_Illegitimate" + file_input_extension, 'r') as f_tlx:
            f_tlx.readline()
            for line in f_tlx:
                if pd.isnull(metadata.loc[metadata['Library'] == label]['MID'].values[0]):
                    len_mid = 0
                else:
                    len_mid = len(metadata.loc[metadata['Library'] == label]['MID'].values[0])
                if pd.isnull(metadata.loc[metadata['Library'] == label]['Primer'].values[0]):
                    len_primer = 0
                else:
                    len_primer = len(metadata.loc[metadata['Library'] == label]['Primer'].values[0])
                distance = int(
                    int(line.split("\t")[9]) - int(line.split("\t")[8])) - (len_mid + len_primer)
                if distance not in distance_dict["illegitimates"][label]:
                    distance_dict["illegitimates"][label][distance] = 1
                else:
                    distance_dict["illegitimates"][label][distance] += 1
                total_junction += 1

        # ADD MAX_DISTANCE AND MAX_DISTANCE_COUNT TO ILLEGITIMATES
        max_distance = 0
        max_distance_count = 0
        if len(distance_dict["illegitimates"][label]) > 0:
            # max_distance
            max_distance = max(distance_dict["illegitimates"][
                               label].iteritems(), key=operator.itemgetter(0))[0]
            # max_distance_count
            max_distance_count = max(distance_dict["illegitimates"][
                                     label].iteritems(), key=operator.itemgetter(1))[1]
        distance_dict["illegitimates"][label]["max_distance"] = max_distance
        distance_dict["illegitimates"][label][
            "max_distance_count"] = max_distance_count

        # ADD MAX_DISTANCE AND MAX_DISTANCE_COUNT TO LEGITIMATES
        max_distance = 0
        max_distance_count = 0
        if len(distance_dict["legitimates"][label]) > 0:
            # max_distance
            max_distance = max(distance_dict["legitimates"][
                               label].iteritems(), key=operator.itemgetter(0))[0]
            # max_distance_count
            max_distance_count = max(distance_dict["legitimates"][
                                     label].iteritems(), key=operator.itemgetter(1))[1]
        distance_dict["legitimates"][label]["max_distance"] = max_distance
        distance_dict["legitimates"][label][
            "max_distance_count"] = max_distance_count
        track_number += 1

    # print(distance_dict)

    max_distance_all = 0
    max_distance_count_all = 0
    # Define max_distance_all and max_distance_count_all
    for legitimate, label_dict in distance_dict.items():
        for key, value in label_dict.items():
            if max_distance_all < int(label_dict[key]['max_distance']):
                max_distance_all = int(label_dict[key]['max_distance'])
            if max_distance_count_all < int(label_dict[key]['max_distance_count']):
                max_distance_count_all = int(
                    label_dict[key]['max_distance_count'])

    # print(max_distance_all)
    # print(max_distance_count_all)

    yrange_count = 0

    for bad_label in global_list:
        label = bad_label.split("/")[-1]

        # MANAGE ARRAYS
        #distance_dict["illegitimates"][label] = manageArray(
        #    distance_dict["illegitimates"][label], max_distance_count_all, method)
        #distance_dict["legitimates"][label] = manageArray(
        #    distance_dict["legitimates"][label], max_distance_count_all, method)

        max_label_leg_ille = max(distance_dict["legitimates"][label]["max_distance_count"],distance_dict["illegitimates"][label]["max_distance_count"])
        distance_dict["illegitimates"][label] = manageArray(
            distance_dict["illegitimates"][label], max_label_leg_ille, method)
        distance_dict["legitimates"][label] = manageArray(
            distance_dict["legitimates"][label], max_label_leg_ille, method)

        #distance_dict["illegitimates"][label] = manageArray(distance_dict["illegitimates"][label], total_junction, method)
        #distance_dict["legitimates"][label] = manageArray(distance_dict["legitimates"][label], total_junction, method)

        distance_dict["illegitimates"][label] = collections.OrderedDict(
            sorted(distance_dict["illegitimates"][label].items()))
        distance_dict["legitimates"][label] = collections.OrderedDict(
            sorted(distance_dict["legitimates"][label].items()))

        # CONCATENATE LEGITIMATE AND ILLEGITIMATE DISTANCE AND COUNTS (USE FOR
        # COLORS)
        xranges_leg = []
        colors_leg = []
        xranges_illeg = []
        colors_illeg = []

        # print(distance_dict["illegitimates"][label])
        for key, value in distance_dict["legitimates"][label].items():
            xranges_leg.append((key, 1))
            colors_leg.append((0.0, 0.0, 0.0, value))

        for key, value in distance_dict["illegitimates"][label].items():
            xranges_illeg.append((key, 1))
            colors_illeg.append((1.0, 0.0, 0.0, value))

        ##############################WORKS##############################
        #xranges=[(29, 1),(30, 1),(31, 1)]
        #yrange=(1.8, 1)
        #colors=[(0.0, 0.0, 0.0, 0.012),(0.0, 0.0, 0.0, 0.008),(0.0, 0.0, 0.0, 0.012)]
        #################################################################

        # print("xranges")
        # print(xranges)
        # print(type(xranges))
        # print("yrange")
        # print(yrange[yrange_count])
        # print(type(yrange[yrange_count]))
        # print("colors")
        # print(colors)
        # print(type(colors))

        # LEGITIMATE DISPLAY
        coll = BrokenBarHCollection(
            xranges_leg, yrange[yrange_count], facecolors=colors_leg, edgecolors=colors_leg)
        ax.add_collection(coll)
        center = yrange[yrange_count][0] + yrange[yrange_count][1] / 2.0
        yticks.append(center)
        yticklabels.append(label + "_legi")
        d[label + "legi"] = xranges_leg

        yrange_count += 1

        # ILLEGITIMATE DISPLAY
        coll = BrokenBarHCollection(xranges_illeg, yrange[
                                    yrange_count], facecolors=colors_illeg, edgecolors=colors_illeg)
        ax.add_collection(coll)
        center = yrange[yrange_count][0] + yrange[yrange_count][1] / 2.0
        yticks.append(center)
        yticklabels.append(label + "_illegi")
        d[label + "illegi"] = xranges_illeg

        yrange_count += 1

    # SET UP DISPLAY BACKGROUND
    ax.axis('tight')
    ax.set_xlim([0, max_distance_all+10])
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    ax.set_xticks(range(0, max_distance_all + 10, 5))
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    plt.xlabel("Position from bait primer")
    plt.axes().xaxis.set_minor_locator(MultipleLocator(1))

    # SET UP TITLE
    plt.title(title_label, fontdict=None, loc='center')
    fig = plt.gcf()
    fig.set_size_inches(22, 10)

    # SAVE DISPLAY
    fig.savefig(dir_post + genome + file_output_extension,
                format='png', bbox_inches='tight')
    # plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])
