# A script allowing to create 4 files to visualize LAM-HTGTS result on a
# proper Circos plot & Karyoplot

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
from Bio import SeqIO
import collections
import copy

# A script allowing to create multiples files usable in input
# visualisation scripts


def usage():
    print('Usage:\n')
    print('\tpython ' + sys.argv[0] + ' -m <metadata file> -g <genome type> -o <output mark> -r <reference fasta file> -p <postprocess directory> -q <legitimate file> [-i <input mark> -s <pool size> -t <percent translocation in circos> -c <percent translocation in karyo> -b <bin size>]')
    print('\t\t-h or --help : display this help')
    print('\t\t-m or --file_metadata : metadata file')
    print('\t\t-g or --genome : only filter librairies results with this genome')
    print('\t\t-o or --output_mark : mark added to the output file')
    print('\t\t-r or --file_reference : reference fasta file')
    print('\t\t-p or --dir_post : postprocess directory')
    print('\t\t-q or --file_legitimate : legitimate positions (bait, prey, start, end, flag)')
    print('\t\t-i or --input_mark : marks from input file')
    print('\t\t-s or --size_pool : the number of bases between two junctions to pool them for illegitimate junctions (Default : 100)')
    print('\t\t-t or --percent_transloc_circos : percentage of translocation insane a size_pool to be display as link on Circos plot (Default : 2.0)')
    print('\t\t-c or --percent_transloc_karyo : percentage of translocation insane a size_pool to be display as link on Karyo plot (Default : 2.0)')
    print('\t\t-b or --bin_size : size of bins to create histogram on Circos plot (Default : 5000000)')


def main(argv):

    pd.options.mode.chained_assignment = None

    file_metadata = ""
    genome = ""
    output_mark = ""
    dir_post = ""
    size_pool = 100
    percent_transloc_circos = 2.0
    percent_transloc_karyo = 2.0
    bin_size = 5000000  # 5Mb
    file_reference = ""
    file_legitimate = ""
    input_mark = ""
    file_input_extension = ""
    file_output_extension = ""

    chromLength_save = {}
    chromLengthStrictBin_save = {}

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'm:g:o:r:p:q:i:s:l:t:c:b:', ['file_metadata=', 'genome=', 'output_mark=',
                                                                        'file_reference=', 'dir_post=', 'file_legitimate=', 'input_mark=', 'size_pool=', 'percent_transloc_circos_circos=', 'percent_transloc_circos_karyo=', 'bin_size=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

##############################OPTIONS##############################

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-m', '--file_metadata'):
            file_metadata = arg
        elif opt in ('-g', '--genome'):
            genome = arg
        elif opt in ('-o', '--output_mark'):
            output_mark = arg
        elif opt in ('-r', '--file_reference'):
            file_reference = arg
        elif opt in ('-p', '--dir_post'):
            dir_post = arg
        elif opt in ('-q', '--file_legitimate'):
            file_legitimate = arg
        elif opt in ('-i', '--input_mark'):
            input_mark = arg
        elif opt in ('-s', '--size_pool'):
            size_pool = arg
        elif opt in ('-t', '--percent_transloc_circos'):
            percent_transloc_circos = arg
        elif opt in ('-c', '--percent_transloc_karyo'):
            percent_transloc_karyo = arg
        elif opt in ('-b', '--bin_size'):
            bin_size = arg
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

    # CHECK MARK
    if output_mark == "":
        print("Error : You have to set a mark !\n")
        usage()
        sys.exit(2)

    # CHECK THE REFERENCE FASTA FILE
    if (file_reference[-3:] != ".fa" and file_reference[-4:] != ".fna" and file_reference[-6:] != ".fasta") or not os.path.exists(file_reference):
        print("Error : The reference fasta file is missing, not .fa, .fna or .fasta !\n")
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

    # CHECK LEGITIME FILE
    if file_legitimate == "" or not os.path.exists(file_legitimate):
        print("Error : You have to set a legitimate file !\n")
        usage()
        sys.exit(2)
    else:
        # READ LEGITIME FILE
        df_legitimate_locus = pd.read_table(
            file_legitimate, sep='\t', header=None)
        for index, row in df_legitimate_locus.iterrows():
            if row[0][0:3] != 'chr':
                print("Error : line." + str(index + 1) +
                      ", col.1 of your legitimate locus file !\n")
                print("Error : Unknown chromosome : " + row[0] + " !\n")
                usage()
                sys.exit(2)
            if row[1][0:3] != 'chr':
                print("Error : line." + str(index + 1) +
                      ", col.2 of your legitimate locus file !\n")
                print("Error : Unknown chromosome : " + row[1] + " !\n")
                usage()
                sys.exit(2)
            try:
                if int(row[2]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.3 of your legitimate locus file !\n")
                    print("Error : Start position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.3 of your legitimate locus file !\n")
                print("Error : Unknown start position : " +
                      str(row[2]) + " !\n")
                usage()
                sys.exit(2)

            try:
                if int(row[3]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.4 of your legitimate locus file !\n")
                    print("Error : End position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.4 of your legitimate locus file !\n")
                print("Error : Unknown end position : " + str(row[3]) + " !\n")
                usage()
                sys.exit(2)

            if int(row[3]) < int(row[2]):
                print("Error : line." + str(index + 1) +
                      " of your legitimate locus file !\n")
                print("Error : End position is smaller than start position !\n")
                usage()
                sys.exit(2)

    # CHECK SIZE POOL
    try:
        size_pool = int(size_pool)
        if size_pool < 1:
            print("Error : Size pool option needs to be more than 0 !\n")
            usage()
            sys.exit(2)
    except:
        print("Error : You have to set an integer to size pool option !\n")
        usage()
        sys.exit(2)

    # CHECK PERCENT TRANSLOCATION CIRCOS
    try:
        percent_transloc_circos = float(percent_transloc_circos)
        if percent_transloc_circos < 0.0:
            print("Error : Percent translocation option needs to be positive !\n")
            usage()
            sys.exit(2)
    except:
        print("Error : You have to set a float to percent translocation option !\n")
        usage()
        sys.exit(2)

    # CHECK PERCENT TRANSLOCATION KARYO
    try:
        percent_transloc_karyo = float(percent_transloc_karyo)
        if percent_transloc_karyo < 0.0:
            print("Error : Percent translocation option needs to be positive !\n")
            usage()
            sys.exit(2)
    except:
        print("Error : You have to set a float to percent translocation option !\n")
        usage()
        sys.exit(2)

    # CHECK BIN SIZE
    try:
        bin_size = int(bin_size)
        if bin_size < 1:
            print("Error : Bin size option needs to be more than 0 !\n")
            usage()
            sys.exit(2)
    except:
        print("Error : You have to set an integer to bin size option !\n")
        usage()
        sys.exit(2)

    # CHECK INPUT MARK HISTORY
    if input_mark == "":
        print("Warning : You will process the raw file !\n")

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
            :-4] + "_" + output_mark + ".csv"
    else:
        file_output_extension = "_" + output_mark + ".csv"

##############################PRINTS##############################

    print('\n-----------------------------------------')
    print('Metadata file : ' + file_metadata)
    print('Genome : ' + genome)
    print('Reference file : ' + file_reference)
    print('Postprocess directory : ' + dir_post)
    print('Legitimate locus file : ' + file_legitimate)
    print('Size pool : ' + str(size_pool))
    print('Percent translocation circos : ' + str(percent_transloc_circos))
    print('Percent translocation karyo : ' + str(percent_transloc_karyo))
    print('Bin size: ' + str(bin_size))
    print('Input file extension: ' + file_input_extension)
    print('Output file extension : ' + file_output_extension)
    print('-----------------------------------------\n')

##############################PROGRAM##############################

    # HIT TABLE INITIALISATION (BINS CREATION FOR HISTOGRAM)
    chromLength_save = getChromLength(file_reference, bin_size)
    
    # HIT TABLE INITIALISATION (BINS CREATION FOR HISTOGRAM)
    chromLengthStrictBin_save = getChromLengthStrictBin(file_reference, size_pool, df_legitimate_locus)

    # LOOP OVER EACH LIBRARIES
    for library in metadata['Library'].tolist():
        print(library)
        chromLength = {}
        chromLength = copy.deepcopy(chromLength_save)
        chromLengthStrictBinPlus = copy.deepcopy(chromLengthStrictBin_save)
        chromLengthStrictBinMinus = copy.deepcopy(chromLengthStrictBin_save)
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

            # MERGE DATAFRAMES
            df = pd.concat([df_legitimate, df_illegitimate], ignore_index=True)
            df = df.sort_values(['Rname', 'Junction'], ascending=[True, True])

            # OUTPUT DISTRIBUTION FILE
            with open(dir_post + library + "/" + library + "_Distribution" + file_output_extension, 'w') as f_link:
                spamwriter = csv.writer(f_link, delimiter='\t')
                spamwriter.writerow(
                    ["Chromosome", "Start", "End", "Strand"])
                for index, row in df.iterrows():
                        spamwriter.writerow(
                            [row['Rname'], row['Junction'], row['Junction'], row['Strand']])

            # FIND NUMBER OF JUNCTION IN BINS ACCORDING TO THE JUNCTION VALUE
            # E.G if bin_size = 100
            # 14950 go to 14900 bin
            # 15500 go to 15500 bin

            i = 0
            total_mutation = 0
            check = False

            # FILL DICTIONNARY FOR CIRCOS HISTOGRAM
            for index, row in df.iterrows():
                if row['Rname'] in chromLength:
                    for key, value in chromLength[row['Rname']].items():
                        # print(row['Junction'])
                        # print(str(key))
                        if int(row['Junction']) <= key and int(row['Junction']) > (key - bin_size):
                            chromLength[row['Rname']][key] += 1
                            total_mutation += 1
                            break
                i += 1
            # TRANSFORM THE NUMBER OF JUNCTION BY PERCENTAGE OF TRANSFORMATION
            for key, value in chromLength.items():
                for key2, value2 in chromLength[key].items():
                    if chromLength[key][key2] != 0:
                        chromLength[key][key2] = float(
                            float(chromLength[key][key2]) / float(total_mutation) * 100)


            # FILL DICTIONNARY FOR KARYO PLOT FREQUENCY
            # ONLY KEEP LEGITIME JUNCTIONS

            for index, row in df_legitimate.iterrows():
                if row['Strand'] == 1:
                    if row['Rname'] in chromLengthStrictBinPlus:
                        for key, value in chromLengthStrictBinPlus[row['Rname']].items():
                            # print(row['Junction'])
                            # print(str(key))
                            if int(row['Junction']) < key and int(row['Junction']) >= (key - size_pool):
                                chromLengthStrictBinPlus[row['Rname']][key] += 1
                                break
                else:
                    if row['Rname'] in chromLengthStrictBinMinus:
                        for key, value in chromLengthStrictBinMinus[row['Rname']].items():
                            # print(row['Junction'])
                            # print(str(key))
                            if int(row['Junction']) < key and int(row['Junction']) >= (key - size_pool):
                                chromLengthStrictBinMinus[row['Rname']][key] += 1
                                break

            # TRANSFORM THE NUMBER OF JUNCTION BY PERCENTAGE OF TRANSFORMATION FOR PLUS
            for key, value in chromLengthStrictBinPlus.items():
                for key2, value2 in chromLengthStrictBinPlus[key].items():
                    if chromLengthStrictBinPlus[key][key2] != 0:
                        chromLengthStrictBinPlus[key][key2] = float(
                            float(chromLengthStrictBinPlus[key][key2]) / float(len(df_legitimate)) * 100)
            # TRANSFORM THE NUMBER OF JUNCTION BY PERCENTAGE OF TRANSFORMATION FOR MINUS
            for key, value in chromLengthStrictBinMinus.items():
                for key2, value2 in chromLengthStrictBinMinus[key].items():
                    if chromLengthStrictBinMinus[key][key2] != 0:
                        chromLengthStrictBinMinus[key][key2] = float(
                            float(chromLengthStrictBinMinus[key][key2]) / float(len(df_legitimate)) * 100)


            # HISTOGRAM FILE CREATION
            Rname_list = []
            Rstart_list = []
            Rend_list = []
            Data_list = []
            dfHisto = pd.DataFrame(columns=['Rname', 'Rstart', 'Rend', 'Data'])
            for key, value in chromLength.items():
                for key_inside, value_inside in value.items():
                    if key_inside - bin_size <= 0:
                        Rname_list.append(key)
                        Rstart_list.append(1)
                        Rend_list.append(key_inside)
                        Data_list.append(value_inside)
                    else:
                        Rname_list.append(key)
                        Rstart_list.append(int(key_inside - (bin_size - 1)))
                        Rend_list.append(key_inside)
                        Data_list.append(value_inside)

            dfHisto = pd.DataFrame({'Rname': pd.Series(Rname_list, dtype=str), 'Rstart': pd.Series(
                Rstart_list, dtype=int), 'Rend': pd.Series(Rend_list, dtype=int), 'Data': pd.Series(Data_list, dtype=float)})

            dfHisto = dfHisto[['Rname', 'Rstart', 'Rend', 'Data']]
            dfHisto = dfHisto.sort_values(
                ['Rname', 'Rstart'], ascending=[True, True])

            # KARYO PLUS FILE CREATION
            Rname_list = []
            Rstart_list = []
            Rend_list = []
            Data_list = []
            Strand_list = []
            dfKaryoFreq = pd.DataFrame(columns=['Rname', 'Rstart', 'Rend', 'Data', 'Strand'])
            for key, value in chromLengthStrictBinPlus.items():
                last_item_key = next(reversed(value))
                for key_inside, value_inside in value.items():
                    if key_inside - size_pool <= 0:
                        Rname_list.append(key)
                        Rstart_list.append(1)
                        Rend_list.append(int(key_inside)-1)
                        Data_list.append(value_inside)
                        Strand_list.append("+")
                    else:
                        if key_inside != last_item_key:
                            Rname_list.append(key)
                            Rstart_list.append(int(key_inside - size_pool))
                            Rend_list.append(int(key_inside)-1)
                            Data_list.append(value_inside)
                            Strand_list.append("+")
                        else:
                            Rname_list.append(key)
                            Rstart_list.append(int(previous_item))
                            Rend_list.append(int(key_inside))
                            Data_list.append(value_inside)
                            Strand_list.append("+")
                        previous_item = key_inside

            for key, value in chromLengthStrictBinMinus.items():
                last_item_key = next(reversed(value))
                for key_inside, value_inside in value.items():
                    if key_inside - size_pool <= 0:
                        Rname_list.append(key)
                        Rstart_list.append(1)
                        Rend_list.append(int(key_inside)-1)
                        Data_list.append(value_inside)
                        Strand_list.append("-")
                    else:
                        if key_inside != last_item_key:
                            Rname_list.append(key)
                            Rstart_list.append(int(key_inside - size_pool))
                            Rend_list.append(int(key_inside)-1)
                            Data_list.append(value_inside)
                            Strand_list.append("-")
                        else:
                            Rname_list.append(key)
                            Rstart_list.append(int(previous_item))
                            Rend_list.append(int(key_inside))
                            Data_list.append(value_inside)
                            Strand_list.append("-")
                        previous_item = key_inside

            dfKaryoFreq = pd.DataFrame({'Rname': pd.Series(Rname_list, dtype=str), 'Rstart': pd.Series(
                Rstart_list, dtype=int), 'Rend': pd.Series(Rend_list, dtype=int), 'Data': pd.Series(Data_list, dtype=float), 'Strand': pd.Series(Strand_list, dtype=str)})

            dfKaryoFreq = dfKaryoFreq[['Rname', 'Rstart', 'Rend', 'Data', 'Strand']]
            dfKaryoFreq = dfKaryoFreq.sort_values(
                ['Strand', 'Rname', 'Rstart'], ascending=[True, True, True])



            # Link modification, keep link above a percentage of translocation
            # PUT df Rname, start and end in array
            # Percentage process
            # Confront each df line against the array
            # DEFAULT : Delete if < 2%
            final_table = []
            i = 0
            j = 0
            nb_translocation = 0
            list_chr = []
            table_position = []
            # This while loop group junctions by bait-prey
            while i < len(df["Qname"]):
                nb_translocation += 1
                if (df["B_Rname"][i] + "-" + df["Rname"][i]) not in list_chr:
                    list_chr.append(df["B_Rname"][i] + "-" + df["Rname"][i])
                    table_position.append(
                        [df["B_Rname"][i] + "-" + df["Rname"][i]])
                table_position[list_chr.index(df["B_Rname"][i] + "-" + df["Rname"][i])].append("-".join([str(df["Junction"][i]), str(
                    df["Rstart"][i]), str(df["Rend"][i]), str(df["Rname"][i]), str(df["B_Rstart"][i]), str(df["B_Rend"][i])]))
                i += 1
            # print(table_position[1])

            # For loop over each bait-prey pair
            table_aux = []
            for i in table_position:
                # print(i)
                # In table_aux append the bait name
                table_aux.append(i[0].split("-")[0])

                # Sort the rest according to the junction value (which is first
                # position)
                for j in i[1:]:
                    # print("*******************************************")
                    # print(table_aux)
                    # print(j)
                    index = 1
                    if len(table_aux[1:]) == 0:
                        #print("INSERT FIRST BECAUSE IT IS EMPTY")
                        table_aux.append(j)
                    else:
                        # Sort by increasing junction value (no need to make a
                        # second sort if junction value are equal, because they
                        # will be group afterwards due to parttern size)
                        if int(j.split("-")[0]) >= int(table_aux[-1].split("-")[0]):
                            #print("INSERT AT THE END")
                            table_aux.append(j)
                        else:
                            #print("INSERT INSIDE")
                            # print(j)
                            for k in table_aux[1:]:
                                if int(j.split("-")[0]) <= int(k.split("-")[0]):
                                    table_aux.insert(index, j)
                                    break
                                index += 1

                dic_size_pool = {}
                nb_transloc = 1
                junction = 0
                min_pos = 0
                max_pos = 0
                chrB = ""
                min_B = 0
                max_B = 0
                j = 1

                # print("**************************************************************************")
                # print(table_aux)
                # print("\n")

                # Create Megajunctions if junction values are close from bin_size value
                # For each chromosome bait
                while j < len(table_aux):
                    # print("-------------------------------------------")
                    # print(j)
                    # print(table_aux[j])
                    # If last
                    if j == len(table_aux) - 1:
                        #print("last cell")
                        # If MegaValues empty
                        if int(min_pos) == 0:
                            #print("NOTHING IN MIN_POS")
                            junction, min_pos, max_pos, chrB, min_B, max_B = firstValues(
                                table_aux[j], junction, min_pos, max_pos, chrB, min_B, max_B)
                        # Update Megavalues
                        else:
                            #print("I got something interesting")
                            junction, min_pos, max_pos, chrB, min_B, max_B = udpateMegajunctionValues(
                                table_aux[j], junction, min_pos, max_pos, chrB, min_B, max_B)
                            #print("CHANGED : "+str(junction)+"-"+str(min_pos)+"-"+str(max_pos)+"-"+chrB+"-"+str(min_B)+"-"+str(max_B))
                        # Write the Megavalue
                        # print("WRITE")
                        # print("-".join([str(junction),str(min_pos),str(max_pos),str(chrB),str(min_B),str(max_B)]))
                        dic_size_pool["-".join([str(abs(junction / nb_transloc)), str(min_pos), str(
                            max_pos), str(chrB), str(min_B), str(max_B)])] = nb_transloc
                        # print("RESET")
                        junction = 0
                        min_pos = 0
                        max_pos = 0
                        chrB = ""
                        min_B = 0
                        max_B = 0
                    # If not last
                    else:
                        #print("not last cell")
                        # If the next one can create a Megavalue
                        # To do so :
                        # |junction-next_junction| <= size_pool

                        if acceptableMergeJunctions(table_aux[j], table_aux[j + 1], int(size_pool)):
                            #print("I got something interesting")
                            junction, min_pos, max_pos, chrB, min_B, max_B = udpateMegajunctionValues(
                                table_aux[j], junction, min_pos, max_pos, chrB, min_B, max_B)
                            #print("CHANGED : "+str(junction)+"-"+str(min_pos)+"-"+str(max_pos)+"-"+chrB+"-"+str(min_B)+"-"+str(max_B))
                            nb_transloc += 1
                        # If the next one can not create a Megavalue
                        else:
                            #print("Cannot catch further")
                            # If MegaValues empty
                            if int(min_pos) == 0:
                                #print("NOTHING IN MIN_POS")
                                junction, min_pos, max_pos, chrB, min_B, max_B = firstValues(
                                    table_aux[j], junction, min_pos, max_pos, chrB, min_B, max_B)
                            # Update Megavalues
                            else:
                                #print("I got something interesting")
                                junction, min_pos, max_pos, chrB, min_B, max_B = udpateMegajunctionValues(
                                    table_aux[j], junction, min_pos, max_pos, chrB, min_B, max_B)
                                #print("CHANGED : "+str(junction)+"-"+str(min_pos)+"-"+str(max_pos)+"-"+chrB+"-"+str(min_B)+"-"+str(max_B))
                            # Write the Megavalue
                            # print("WRITE")
                            #print("NB_TRANSLOC : "+str(nb_transloc))
                            # print("-".join([str(junction),str(min_pos),str(max_pos),str(chrB),str(min_B),str(max_B)]))
                            dic_size_pool["-".join([str(abs(junction / nb_transloc)), str(min_pos), str(
                                max_pos), str(chrB), str(min_B), str(max_B)])] = nb_transloc
                            junction = 0
                            min_pos = 0
                            max_pos = 0
                            chrB = ""
                            min_B = 0
                            max_B = 0
                            # print("RESET")

                            nb_transloc = 1

                    # print(dic_size_pool)
                    j += 1

                # print("----------------CHECK-------------------\n")
                # print(nb_translocation)
                # print(dic_size_pool)

                for key, value in dic_size_pool.items():
                    dic_size_pool[key] = round(
                        float(int(value) * 100 / float(nb_translocation)), 2)

                del table_aux[1:]
                table_aux.append(dic_size_pool)
                final_table.append(table_aux)
                index = 0
                table_aux = []

            # OUTPUT

            # i=0
            # while i<len(df["Qname"]):
            #    df["Qname"][i]=df["Qname"][i].split(":")[0]+":"+df["Qname"][i].split(":")[-2]+":"+df["Qname"][i].split(":")[-1]
            #    i+=1
            # dfNormal=df
            #dfNormal[["Rname", "Rstart", "Rend", "Qname"]].to_csv(output+"normal.csv", sep='\t', encoding='utf-8', index=False)

            dfHisto[["Rname", "Rstart", "Rend", "Data"]].to_csv(
                dir_post + library + "/" + library + "_Histo" + file_output_extension, sep='\t', encoding='utf-8', index=False)

            dfKaryoFreq[["Rname", "Rstart", "Rend", "Data", "Strand"]].to_csv(
                dir_post + library + "/" + library + "_Karyoplot_freq" + file_output_extension, sep='\t', encoding='utf-8', index=False)

            chrP = ""
            with open(dir_post + library + "/" + library + "_Link" + file_output_extension, 'w') as f_link:
                spamwriter = csv.writer(f_link, delimiter='\t')
                spamwriter.writerow(
                    ["Rname", "Rstart", "Rend", "B_Rname", "B_Rstart", "B_Rend"])
                for j in final_table:
                    for key, value in j[1].items():
                        if float(value) > float(percent_transloc_circos):
                            chrP = key.split("-")[3]
                            min_pos = key.split("-")[1]
                            max_pos = key.split("-")[2]
                            chrB = j[0]
                            min_B = key.split("-")[4]
                            max_B = key.split("-")[5]
                            spamwriter.writerow(
                                [chrP, min_pos, max_pos, chrB, min_B, max_B])

            with open(dir_post + library + "/" + library + "_Karyoplot_link" + file_output_extension, 'w') as f_link:
                spamwriter = csv.writer(f_link, delimiter='\t')
                spamwriter.writerow(
                    ["Rname", "Junction", "Rstart", "Rend", "B_Rname", "B_Rstart", "B_Rend", "value"])
                for j in final_table:
                    for key, value in j[1].items():
                        if float(value) > float(percent_transloc_karyo):
                            chrP = key.split("-")[3]
                            junction = key.split("-")[0]
                            min_pos = key.split("-")[1]
                            max_pos = key.split("-")[2]
                            chrB = j[0]
                            min_B = key.split("-")[4]
                            max_B = key.split("-")[5]
                            spamwriter.writerow(
                                [chrP, junction, min_pos, max_pos, chrB, min_B, max_B, value])
        # sys.exit()

##############################FUNCTIONS##############################

# ENABLE ITERATION ON DATAFRAME


def iterate(iterable):
    iterator = iter(iterable)
    item = iterator.next()

    for next_item in iterator:
        yield item, next_item
        item = next_item

    yield item, None

# Initialize bin for histogram representation in circos plot


def getChromLength(file_reference, bin_size):
    chromLength = {}
    chrom_dict = {}
    chrom_dict = {
        record.id: len(record.seq)
        for record in SeqIO.parse(file_reference, 'fasta')
        if record.id not in chrom_dict
    }
    # Bin creation by chromosome starting from the end of each chromosome,
    # create bin_size bin size
    for key, value in chrom_dict.items():
        max_value = value
        chromLength[key] = {}
        while int(max_value) > bin_size:
            chromLength[key][int(max_value)] = 0
            max_value -= bin_size
        chromLength[key][int(max_value)] = 0
        chromLength[key] = collections.OrderedDict(
            sorted(chromLength[key].items()))
    return chromLength

# Initialize bin for frequency representation in karyo plot

def getChromLengthStrictBin(file_reference, bin_size, df_legitimate_locus):
    chromLength = {}
    chrom_dict = {}
    for index, row in df_legitimate_locus.iterrows():
        key = row[1]
        start = row[2]
        end = row[3]
        chromLength[key] = {}
        while start+bin_size < end:
            start += bin_size
            chromLength[key][start] = 0
        chromLength[key][start] = 0
        chromLength[key][end] = 0
        chromLength[key] = collections.OrderedDict(
            sorted(chromLength[key].items()))
    return chromLength

# Extend Megajunction values to get the hugest junction (merge junctions
# start and end positions)


def udpateMegajunctionValues(table, junction, min_pos, max_pos, chrB, min_B, max_B):
        # if int(table.split("-")[0])>int(junction) or int(junction)==0:
        #print("junction changed")
        #print("junction : "+str(junction))
        #print("conflict : "+str(table.split("-")[0]))
        # junction=int(table.split("-")[0])
    junction += int(table.split("-")[0])
    if int(table.split("-")[1]) < int(min_pos) or int(min_pos) == 0:
        #print("min_pos changed")
        #print("min_pos : "+str(min_pos))
        #print("conflict : "+str(table.split("-")[1]))
        min_pos = int(table.split("-")[1])
    if int(table.split("-")[2]) > int(max_pos):
        #print("max_pos changed")
        #print("max_pos : "+str(max_pos))
        #print("conflict : "+str(table.split("-")[2]))
        max_pos = int(table.split("-")[2])
    if int(table.split("-")[4]) < int(min_B) or int(min_B) == 0:
        #print("min_B changed")
        #print("min_B : "+str(min_B))
        #print("conflict : "+str(table.split("-")[4]))
        min_B = int(table.split("-")[4])
        if str(table.split("-")[3]) != str(chrB) and str(chrB) != "":
            #print("On va avoir un probleme chef !")
            sys.exit()
        else:
            chrB = str(table.split("-")[3])
    if int(table.split("-")[5]) > int(max_B):
        #print("max_B changed")
        #print("max_B : "+str(max_B))
        #print("conflict : "+str(table.split("-")[5]))
        max_B = int(table.split("-")[5])
        if str(table.split("-")[3]) != str(chrB) and str(chrB) != "":
            #print("On va avoir un probleme chef !")
            sys.exit()
        else:
            chrB = str(table.split("-")[3])
    return junction, min_pos, max_pos, chrB, min_B, max_B

# Initialize Megajunction values if none

def firstValues(table, junction, min_pos, max_pos, chrB, min_B, max_B):
    junction = int(table.split("-")[0])
    min_pos = int(table.split("-")[1])
    max_pos = int(table.split("-")[2])
    chrB = str(table.split("-")[3])
    min_B = int(table.split("-")[4])
    max_B = int(table.split("-")[5])
    #print("junction : "+str(junction))
    #print("min_pos : "+str(min_pos))
    #print("max_pos : "+str(max_pos))
    #print("chrB : "+str(chrB))
    #print("min_B : "+str(min_B))
    #print("max_B : "+str(max_B))
    return junction, min_pos, max_pos, chrB, min_B, max_B

# Validate or not the Megajunction possibility


def acceptableMergeJunctions(actual, next, size_pool):
    check_size_pool = False
    check_bait = False
    check_prey = False

    #print("actual : "+actual)
    #print("next   : " +next)
    # Not necessary to check with absolute value because it is sorted

    # Check for junction size_pool
    if int(next.split("-")[0]) - int(actual.split("-")[0]) <= size_pool:
        #print("size_pool checked")
        check_size_pool = True

    # check for bait size_pool (start_end compatibility)
    if int(actual.split("-")[2]) >= int(next.split("-")[2]) + size_pool:
        if int(actual.split("-")[1]) <= int(next.split("-")[2]) + size_pool:
            #print("bait checked")
            check_bait = True
    else:
        if int(actual.split("-")[2]) >= int(next.split("-")[1]) - size_pool:
            #print("bait checked")
            check_bait = True

    # check for prey size_pool (start_end compatibility) (not really necessary
    # because of size_pool_check but still here)
    if int(actual.split("-")[5]) >= int(next.split("-")[5]) + size_pool:
        if int(actual.split("-")[4]) <= int(next.split("-")[5]) + size_pool:
            #print("prey checked")
            check_prey = True
    else:
        if int(actual.split("-")[5]) >= int(next.split("-")[4]) - size_pool:
            #print("prey checked")
            check_prey = True

    # All check completed
    if check_size_pool and check_bait and check_prey:
        return True
    else:
        return False


if __name__ == '__main__':
    main(sys.argv[1:])
