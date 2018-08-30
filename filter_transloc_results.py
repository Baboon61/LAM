# A script allowing to filter TranslocPipeline results calling the
# TranslocFilter.pl function on the metadata.txt file

##############################IMPORTS##############################

import os
import sys
import getopt
import pandas as pd

# A script allowing to filter TranslocPipeline results calling the
# TranslocFilter.pl function on the metadata.txt file

# All options are coming from the https://robinmeyers.github.io/transloc_pipeline/thedocs.html
# But no documentation on how to use these options for now
# Here are the options unaligned, baitonly, misprimed, freqcut, largegap, mapqual, breaksite, sequential, repeatseq, duplicate
# Default : f.unaligned=1 f.baitonly=1 f.misprimed=L10 f.freqcut=1 f.largegap=G30 f.mapqual=1 f.breaksite=1 f.sequential=1 f.repeatseq=1 f.duplicate=1
# Two types of option, binary option (1 or 0) and value(L -> lower than and G -> greater than, just the value if you want to filter one single value)
# These script will works on all your metadata libraries results


def usage():
    print('Usage:')
    print('\tpython ' + sys.argv[0] + ' -f <metadata file> -t <results directory> -o <output mark> -g <genome type> [-i <input mark> --unaligned <value> --baitonly <value> --uncut <value> --misprimed <value> --freqcut <value> --largegap <value> --mapqual <value> --breaksite <value> --sequential <value> --repeatseq <value> --duplicate <value>]')
    print('\t\t-h or --help : display this help')
    print('\t\t-m or --file_metadata : metadata file')
    print('\t\t-t or --dir_results : results directory')
    print('\t\t-o or --output_mark : mark added to the output file name')
    print('\t\t-g or --genome : only filter libraries results with this genome')
    print('\t\t-i or --input_mark : marks from input file')
    print('\t\t--unaligned : No OCS alignments (0 or 1)')
    print('\t\t--baitonly : Bait alignment is either the only alignent in the OCS or only followed by adapter alignment (0 or 1)')
    print('\t\t--uncut : Bait alignment runs greater than some number of bases past the cutsite (L or G, then your value)')
    print('\t\t--misprimed : Bait alignment runs fewer than some number of bases past the primer (L or G, then your value)')
    print('\t\t--freqcut : Restriction enzyme site within some number of bases from the junction (0 or 1)')
    print('\t\t--largegap : More than some number of bases between the bait and prey alignments (L or G, then your value)')
    print('\t\t--mapqual : OCS had a competing prey junction (0 or 1)')
    print('\t\t--breaksite : Prey alignment maps into non-endogenous breaksite cassette (0 or 1)')
    print('\t\t--sequential : Junction occurs downstream on read from first bait-prey junction (0 or 1)')
    print('\t\t--repeatseq : Repeated sequences (0 or 1)')
    print('\t\t--duplicate : Duplicate junctions (0 or 1)')
    print("INFO : To not filter just don't put it in options")


def main(argv):

    command_line = ""
    options_line = ""
    file_metadata = ""
    dir_results = ""
    output_mark = ""
    genome = ""
    input_mark = ""
    file_input_extension = ""
    file_output_extension = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'm:t:o:g:i:', ['file_metadata=', 'dir_results=', 'output_mark=', 'genome=', 'input_mark=', 'unaligned=',
                                                                'baitonly=', 'uncut=', 'misprimed=', 'freqcut=', 'largegap=', 'mapqual=', 'breaksite=', 'sequential=', 'repeatseq=', 'duplicate=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

##############################OPTIONS##############################

    if not opts:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-m', '--file_metadata'):
            file_metadata = arg
        elif opt in ('-t', '--dir_results'):
            directory_path = arg
        elif opt in ('-o', '--output_mark'):
            output_mark = arg
        elif opt in ('-g', '--genome'):
            genome = arg
        elif opt in ('-i', '--input_mark'):
            input_mark = arg
        elif opt in ('--unaligned'):
            options_line += " f.unaligned=" + arg
        elif opt in ('--baitonly'):
            options_line += " f.baitonly=" + arg
        elif opt in ('--uncut'):
            options_line += " f.uncut=" + arg
        elif opt in ('--misprimed'):
            options_line += " f.misprimed=" + arg
        elif opt in ('--freqcut'):
            options_line += " f.freqcut=" + arg
        elif opt in ('--largegap'):
            options_line += " f.largegap=" + arg
        elif opt in ('--mapqual'):
            options_line += " f.mapqual=" + arg
        elif opt in ('--breaksite'):
            options_line += " f.breaksite=" + arg
        elif opt in ('--sequential'):
            options_line += " f.sequential=" + arg
        elif opt in ('--repeatseq'):
            options_line += " f.repeatseq=" + arg
        elif opt in ('--duplicate'):
            options_line += " f.duplicate=" + arg
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

    # CHECK RESULTS DIRECTORY
    if not os.path.exists(dir_results):
        print("Error : You have to set a results directory !\n")
        usage()
        sys.exit(2)
    else:
        if dir_results[-1] != "/":
            dir_results += "/"

    # FILTER METADATA FILE IF GENOME INPUT
    if genome != "":
        metadata = metadata.loc[metadata['Assembly'] == genome]
        if metadata.empty:
            print("Error : This assembly does not exist in the metadata file !\n")
            usage()
            sys.exit(2)

    # CHECK INPUT MARKS HISTORY
    if input_mark == "":
        print("Warning : You will filter the raw file !\n")

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
    for i in metadata['Library'].tolist():
        if os.path.exists(dir_results + i + "/" + i + file_input_extension):
            check_input_mark = True
    if not check_input_mark:
        print("Error : Your input marks can not localize a good input file !\n")
        usage()
        sys.exit(2)

    # SELECT OUTPUT FILES
    if file_input_extension != "":
        file_output_extension = file_input_extension[
            :-4] + "_" + output_mark + ".tlx"
    else:
        file_output_extension = "_" + output_mark + ".tlx"

##############################PRINTS##############################

    print('\n-----------------------------------------')
    print('Metadata file : ' + file_metadata)
    print('Results directory : ' + dir_results)
    print('Genome : ' + genome)
    print('Input file extension: ' + file_input_extension)
    print('Output file extension : ' + file_output_extension)
    print('-----------------------------------------\n')

##############################PROGRAM##############################

    # LOOP OVER EACH LIBRARIES
    for i in metadata['Library'].tolist():
        # CHECK DIRECTORY EXISTS
        if not os.path.exists(dir_results + i):
            print("Warning : " + dir_results +
                  " does not contains {" + i + "}")
            print("Warning :  {" + i + "} will not be filtered")
        else:
            # CHECK INPUT FILE EXISTS
            if os.path.exists(directory_path + i + "/" + i + file_input_extension):
                # CREATE THE COMMAND LINE
                command_line = "TranslocFilter.pl " + dir_results + i + "/" + i + file_input_extension + " " + \
                    dir_results + i + "/" + i + file_output_extension + \
                    " --filters " + "'" + options_line[1:] + "'"
                print(command_line)
                os.system(command_line)
            else:
                print("Warning : " + dir_results + i + "/" + i +
                      file_input_extension + " does not exist")
                print("Warning :  {" + dir_results + i + "/" + i +
                      file_input_extension + "} will not be filtered")

if __name__ == '__main__':
    main(sys.argv[1:])
