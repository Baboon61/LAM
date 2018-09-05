# A script allowing to run CSReport software on all the librairies of
# metadata.txt

##############################IMPORTS##############################

import os
import sys
import shutil
import getopt
import pandas as pd


def usage():
    print('Usage:')
    print('\tpython ' + sys.argv[0] + ' -m <metadata file> -g <genome type> -r <reference fasta file> -p <postprocess directory> -f <annotation file> -a <from annotation> -b <to annotation> [-i <mark input> -u <unlink> -s <cluster size> -z <reverse CSReport>]')
    print('\t\t-h or --help : display this help')
    print('\t\t-m or --file_metadata : metadata file')
    print('\t\t-g or --genome : only filter librairies results with this genome')
    print('\t\t-r or --file_reference : reference fasta file')
    print('\t\t-p or --dir_post : postprocess directory')
    print('\t\t-f or --file_annot : annotation bed file (locus, start, end)')
    print("\t\t-a or --from_annot : search junction from locis, separated by comma, write 'all' to confront from all locis")
    print("\t\t-b or --to_annot : search junction to locis, separated by comma, write 'all' to confront to all locis")
    print('\t\t-i or --input_mark : mark from input file')
    print('\t\t-u or --unlink : unlink the from_annot1/to_annot1, from_annot2/to_annot2 process and allows to search for from_annot1/to_annot2 junctions')
    print('\t\t-s or --cluster_size : count number of clustered junctions (Default : 2)')
    print('\t\t-z or --reverse_CSR : do normal + reverse CSReport (change from_annot and to_annot)')


def main(argv):

    path_CSReport = os.path.dirname(os.path.abspath(__file__)) + "/CSReport/"
    file_metadata = ""
    genome = ""
    dir_post = ""
    file_annot = ""
    from_annot = ""
    to_annot = ""
    unlink = False
    file_reference = ""
    cluster_size = 2
    input_mark = ""
    reverse_CSR = False
    file_input_extension = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'm:g:r:p:f:a:b:i:us:z', ['file_metadata=', 'genome=', 'file_reference=', 'dir_post=',
                                                                          'file_annot=', 'from_annot=', 'to_annot=', 'input_mark=', 'unlink', 'cluster_size=', 'reverse_CSR=', 'help'])
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
        elif opt in ('-g', '--genome'):
            genome = arg
        elif opt in ('-r', '--file_reference'):
            file_reference = arg
        elif opt in ('-p', '--dir_post'):
            dir_post = arg
        elif opt in ('-f', '--file_annot'):
            file_annot = arg
        elif opt in ('-a', '--from_annot'):
            from_annot = arg
        elif opt in ('-b', '--to_annot'):
            to_annot = arg
        elif opt in ('-i', '--input_mark'):
            input_mark = arg
        elif opt in ('-u', '--unlink'):
            unlink = True
        elif opt in ('-s', '--cluster_size'):
            cluster_size = arg
        elif opt in ('-z', '--reverse_CSR'):
            reverse_CSR = True
        else:
            print("Error : Bad option")
            usage()
            sys.exit(2)

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

    # CHECK REFERENCE FILE
    if (file_reference[-3:] != ".fa" and file_reference[-4:] != ".fna" and file_reference[-6:] != ".fasta") or not os.path.exists(file_reference):
        print("Error : The reference file is not fasta or does not exist !\n")
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

    # CHECK ANNOTATION FILE
    if file_annot == "" or not os.path.exists(file_annot):
        print("Error : You have to set an annotation file !\n")
        usage()
        sys.exit(2)
    else:
        # READ ANNOTATION FILE
        df_annot = pd.read_table(file_annot, sep='\t', header=None)
        for index, row in df_annot.iterrows():
            try:
                if int(row[1]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.2 of your annotation file !\n")
                    print("Error : Start position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.2 of your annotation file !\n")
                print("Error : Unknown start position : " +
                      str(row[1]) + " !\n")
                usage()
                sys.exit(2)
            try:
                if int(row[2]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.3 of your annotation file !\n")
                    print("Error : End position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.3 of your annotation file !\n")
                print("Error : Unknown end position : " + str(row[2]) + " !\n")
                usage()
                sys.exit(2)
            if int(row[2]) < int(row[1]):
                print("Error : line." + str(index + 1) +
                      " of your annotation file !\n")
                print("Error : End position is smaller than start position !\n")
                usage()
                sys.exit(2)

    # CHECK FROM ANNOTATION OPTION
    if from_annot != "":
        if 'all' in from_annot:
            list_from_annot = df_annot[0].tolist()
        else:
            list_from_annot = from_annot.split(",")
            for i in list_from_annot:
                if i not in df_annot[0].tolist():
                    print(
                        "Error : {" + i + "} does not exist in the annotation file !\n")
                    usage()
                    sys.exit(2)

    # CHECK TO ANNOTATION OPTION
    if to_annot != "":
        if 'all' in to_annot:
            list_to_annot = df_annot[0].tolist()
        else:
            list_to_annot = to_annot.split(",")
            for i in list_to_annot:
                if i not in df_annot[0].tolist():
                    print(
                        "Error : {" + i + "} does not exist in the annotation file !\n")
                    usage()
                    sys.exit(2)

    # CHECK UNLINK
    if not unlink:
        if 'all' in from_annot or 'all' in to_annot:
            print("Warning : Default unlink if -a all or -b all !")
            unlink = True
        else:
            if len(list_from_annot) != len(list_to_annot):
                print(
                    "Error : The number of element in from_annotation and to_annotation have to be the same if you link them !")
                usage()
                sys.exit(2)
            unlink = False

    # CHECK CLUSTER SIZE
    try:
        cluster_size = int(cluster_size)
        if cluster_size < 1:
            print("Error : Cluster size option needs to be more than 0 !\n")
            usage()
            sys.exit(2)
    except:
        print("Error : You have to set an integer to cluster size option !\n")
        usage()
        sys.exit(2)

    # CHECK INPUT MARKS HISTORY
    if input_mark == "":
        print("Warning : You will process the raw file !\n")

    # SELECT INPUT FILES
    if input_mark == "":
        file_input_extension = ".tlx"
    else:
        file_input_extension = "_" + "_".join(input_mark.split(",")) + ".fasta"

##############################PRINTS##############################

    print('\n-----------------------------------------')
    print('Metadata file : ' + file_metadata)
    print('Genome : ' + genome)
    print('Reference fasta file : ' + file_reference)
    print('Postprocess directory : ' + dir_post)
    print('Annotation file : ' + file_annot)
    print('From annotation : ' + ','.join(list_from_annot))
    print('To annotation : ' + ','.join(list_to_annot))
    if unlink:
        print('Unlink : Yes')
    else:
        print('Unlink : No')
    print('Cluster size : ' + str(cluster_size))
    if reverse_CSR:
        print('Reverse CSReport : Yes')
    else:
        print('Reverse CSReport : No')
    print('Input file extension: ' + file_input_extension)
    print('-----------------------------------------\n')

##############################PROGRAM##############################

    # IMPORT CSReport
    from CSReport import CSReport_run, CSReport_summary, histogramStructures, sketchBP, evaluateDiversity, summarizeMotifs

    # LOOP OVER EACH LIBRARIES
    for library in metadata['Library'].tolist():
        if not os.path.exists(dir_post + library):
            print("Warning : " + dir_post +
                  " does not contains {" + library + "}")
            print("Warning :  {" + library + "} will not be filtered")
        else:
            # CHECKK CSReport FOLDER
            os.mkdir(dir_post + library + "/CSReport" +
                     file_input_extension[0:-6])
            for j in list_from_annot:
                for k in list_to_annot:
                    print(library, j, k)
                    path_to_fasta = dir_post + library + "/"
                    fasta = library + "_" + \
                        k.replace("/", "-").replace(" ",
                                                    "+").replace("'", "!") + file_input_extension
                    # CHECK SEQUENCES FILE
                    if os.stat(path_to_fasta + fasta).st_size != 0:
                        if not os.path.exists(path_to_fasta + fasta):
                            print("Warning : " + path_to_fasta +
                                  fasta + " does not exist")
                            print(
                                "Warning :  {" + path_to_fasta + fasta + "} will not be reported")
                        else:
                            print("NORMAL")
                            CSReport_run(Seq=path_to_fasta + fasta[0:-6], Ref=file_reference[
                                         0:-6], Region1=j, Region2=k, Form='fasta')
                            if os.path.exists(path_to_fasta + fasta[0:-6] + "/" + fasta[0:-6] + "_readsJunction.csv"):
                                df_junction = pd.read_table(
                                    path_to_fasta + fasta[0:-6] + "/" + fasta[0:-6] + "_readsJunction.csv", sep="\t", header=0)
                                # CREATE THE COMMAND LINE
                                if len(df_junction.index) > 0:
                                    current_path = os.getcwd()
                                    os.chdir(path_to_fasta)

                                    # FOR NOW DISABLE OTHER FUNCTIONS

                                    # CSReport_summary(fasta[0:-6],cluster=cluster_size)
                                    # Structures=histogramStructures(Seq=fasta[0:-6],bw=False)
                                    # sketchBP(Seq=fasta[0:-6],Ref=file_reference[0:-6],Region1=j,Region2=k,density2D=False)

                                    # evaluateDiversity(Seq=fasta[0:-6])
                                    # summarizeMotifs(Seq=fasta[0:-6],Ref=file_reference[0:-6],Region1=j,Region2=k)

                                    os.chdir(current_path)
                                    shutil.move(path_to_fasta + fasta[0:-6], path_to_fasta + "CSReport" + file_input_extension[0:-6] + "/" + j.replace(
                                        "/", "-").replace(" ", "+").replace("'", "!") + "_" + k.replace("/", "-").replace(" ", "+").replace("'", "!") + "_" + fasta[0:-6])
                                else:
                                    shutil.rmtree(path_to_fasta + fasta[0:-6])

                            # IF REVERSE CSReport IS SET
                            if reverse_CSR and j != k:
                                print("REVERSE")
                                CSReport_run(Seq=path_to_fasta + fasta[0:-6], Ref=file_reference[
                                             0:-6], Region1=k, Region2=j, Form='fasta')
                                if os.path.exists(path_to_fasta + fasta[0:-6] + "/" + fasta[0:-6] + "_readsJunction.csv"):
                                    df_junction = pd.read_table(
                                        path_to_fasta + fasta[0:-6] + "/" + fasta[0:-6] + "_readsJunction.csv", sep="\t", header=0)
                                    if len(df_junction.index) > 0:
                                        current_path = os.getcwd()
                                        os.chdir(path_to_fasta)

                                        # FOR NOW DISABLE OTHER FUNCTIONS

                                        # CSReport_summary(fasta[0:-6],cluster=cluster_size)
                                        # Structures=histogramStructures(Seq=fasta[0:-6],bw=False)
                                        # sketchBP(Seq=fasta[0:-6],Ref=file_reference[0:-6],Region1=k,Region2=j,density2D=False)

                                        # evaluateDiversity(Seq=fasta[0:-6])
                                        # summarizeMotifs(Seq=fasta[0:-6],Ref=file_reference[0:-6],Region1=k,Region2=j)

                                        os.chdir(current_path)
                                        shutil.move(path_to_fasta + fasta[0:-6], path_to_fasta + "CSReport" + file_input_extension[0:-6] + "/" + k.replace(
                                            "/", "-").replace(" ", "+").replace("'", "!") + "_" + j.replace("/", "-").replace(" ", "+").replace("'", "!") + "_" + fasta[0:-6])
                                    else:
                                        shutil.rmtree(
                                            path_to_fasta + fasta[0:-6])

            # MOVE ALL SEQUENCE FILES
            print("mv " + path_to_fasta + "*" + file_input_extension + " " +
                  path_to_fasta + "CSReport" + file_input_extension[0:-6] + "/")
            os.system("mv " + path_to_fasta + "*" + file_input_extension + " " +
                      path_to_fasta + "CSReport" + file_input_extension[0:-6] + "/")

if __name__ == '__main__':
    main(sys.argv[1:])
