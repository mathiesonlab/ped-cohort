"""
Main file to find minimum pedigree from an IBD cohort to a source
data.
Authors: Alton Wiggers
Date: 6/7/21
"""

#python imports
import argparse
import sys
import time
import signal
import os #used for testing

#local imports
import IBD
from PedigreeTree import PedigreeTree
from AncestorNode import AncestorNode

default_timeout = 0

class Parser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('\nerror: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args(description): #argument parsing

    parser = Parser()
    parser.add_argument("struct_filename", \
        help="input .txt file with pedigree structure")
    parser.add_argument("germ_filename", \
        help="input GERMLINE .match file")
    parser.add_argument("-o", "--output_filename", \
        help="output .txt file with minimum pedigree structure")
    parser.add_argument("-p", "--pedigree_filenames", nargs=2, \
        help="input and output file names for an associated .ped file (2 arg format: -p [input_ped] [output_ped]")
    parser.add_argument("-i", "--ibd", nargs=4, \
        help="a specific ibd to choose (4 arg format: -i [chr] [start] [end] [indv])")
    parser.add_argument("-s", "--source", \
        help="a specific source to choose (recommended to use with -i). Please use + instead of & for couples.")
    parser.add_argument("-t", "--timeout", type=int, \
        help="set a number of seconds allowed to search for paths for an individual source")
    parser.add_argument("-q", "--quiet", action="store_true", \
        help="supress terminal output")

    mandatories = ["germ_filename", "struct_filename"]


    args = parser.parse_args()

    return args

class TimeoutException(Exception): #timeout exception class
    pass

def timeout_handler(signum, frame): #timeout signal handler
    raise TimeoutException


def main():

    #parse arguments
    args = parse_args("pedigree args")

    #set timeout
    timeout=default_timeout
    if args.timeout != None:
        timeout = args.timeout

    # construct pedigree data structure
    ped = PedigreeTree(args.struct_filename)
    left_out = []

    # construct IBDs data structures (type: list)
    IBDs = IBD.get_IBDs(args.germ_filename, left_out) # toggle to leave out

    # assign IBDs to individuals in pedigree
    IBD.ibd_to_indvs(IBDs, ped)

    #testing code
    # python3 ped-cohort.py input/extended_pedigree_final.fam output/trial1_chr21_amr_geno_germline.match -o output/trial1_chr21_cohort.fam -p output/trial1_chr21_amr_geno.ped output/trial1_chr21_amr_cohort.ped -t 10
    '''
    i = 0
    #info_file = open("recon_info.txt","w")
    #info_file.close()
    for selected_ibd in IBDs:
        starting_indvs = []
        for indv in selected_ibd.get_indvs():
            starting_indvs.append(indv)
        list_options = find_min_pedigree(ped,starting_indvs,args.source,timeout,args.quiet)
        
        if args.output_filename != None:
            for option in list_options:
                write_to_file(args.output_filename,ped,option[1],args.quiet)
                status = os.system("Rscript visualizePed.R " + str(args.output_filename))
                if status != 0:
                    print("failed early with status: " + str(status))
                    exit()
        
        if args.pedigree_filenames != None:
            info_file = open("recon_info.txt","a")
            info_file.write(selected_ibd.id + "\n")
            info_file.close()
            for option in list_options:
                write_to_file(args.output_filename,ped,option[1],args.quiet)
                create_ped_file(args.pedigree_filenames[0], args.pedigree_filenames[1], option[1],args.quiet)
                os.system("germline -input output/trial1_chr21_amr_cohort.ped output/trial1_chr21_amr.map -haploid -output output/trial1_chr21_amr_cohort_germline")
                os.system("python3 match2json.py -g output/trial1_chr21_amr_cohort_germline.match -s output/trial1_chr21_cohort.fam  -m output/trial1_chr21_amr.map -p output/trial1_chr21_amr_cohort.ped -j output/trial1_chr21_amr_cohort.json")
                os.system("PYTHONHASHSEED=1833 python3 thread.py -g output/trial1_chr21_amr_cohort_germline.match -s output/trial1_chr21_cohort.fam -m output/trial1_chr21_amr.map -j output/trial1_chr21_amr_cohort.json -p output/trial1_chr21_amr_recon.ped")
                ped_file = open("output/trial1_chr21_amr_recon.ped","r")
                didReconstruct = False
                for line in ped_file:
                    if line.strip():
                        words = line.split()
                        for i in range(6,len(words)):
                            if words[i] != '?':
                                didReconstruct = True
                                break
                        if didReconstruct:
                            break
                ped_file.close()
                info_file = open("recon_info.txt","a")
                info_file.write(option[0] + " "  + str(starting_indvs) + " " + str(len(option[1])) + " "+ str(didReconstruct) +"\n")
                info_file.close()
        print("completed IBD " + str(i))
        i+=1
    
    print("all IBDs checked successfully")
    exit()
    '''

    
    selected_ibd = None
    if args.ibd == None: #check for ibd selection at command line
        print("chromosome start end indv.haplotype")
        print("-----------------------------------")
        for i in range(len(IBDs)): # print options
            print("[" + str(i) + "] " + IBDs[i].id + " (cohort size: " + str(len(IBDs[i].get_indvs())) + ")")
        
        while selected_ibd == None: #continue prompt until a valid input is received
            user_in = input("Please select an IBD from the list above: ")
            # TODO make multiple selections?
            if not user_in.isdigit() or int(user_in) < 0 or int(user_in) >= len(IBDs):
                print("Invalid input. Please input a number from the lists above.")
            else:
                selected_ibd = IBDs[int(user_in)]
    else: #using preselected ibd
        ibd_id = args.ibd[0] + " " + args.ibd[1] + " " + args.ibd[2] + " " + args.ibd[3]
        for i in range(len(IBDs)):
            if IBDs[i].id == ibd_id:
                selected_ibd = IBDs[i]
                break
        if selected_ibd == None:
            print("Could not find specified IBD")
            exit()
    
    
    """
    ibd_set = {}
    for ibd in IBDs:
        id_elems = ibd.id.split()
        id = id_elems[0] + " " + id_elems[1] + " " + id_elems[2]
        if id in ibd_set:
            ibd_set[id] = ibd_set[id] + [id_elems[3]]
        else:
            ibd_set[id] = [id_elems[3]]
    for id in ibd_set:
        if len(ibd_set[id]) > 1:
            print(id + " has multiple instances: ")
            for mem in ibd_set[id]:
                full_id = id+" "+mem
                for ibd in IBDs:
                    if ibd.id == full_id:
                        print("\t"+str(ibd.get_indvs()))
    """

    starting_indvs = []
    output_list = []

    #find minimum pedigree for all relevant sources
    for indv in selected_ibd.get_indvs():
        starting_indvs.append(indv)
    list_options = find_min_pedigree(ped,starting_indvs,args.source,timeout,args.quiet)
    
    selection_numbers = []

    if len(list_options) == 0:
        print("no pedigree options were found")
        exit()
    elif args.source == None: #check for source selection at command line
        for i in range(len(list_options)): #print options
            print("[" + str(i) + "] source: " + list_options[i][0] + " pedigree size: " + str(len(list_options[i][1])))
    else: #using preselected source
        print("source: " + list_options[0][0] + " pedigree size: " + str(len(list_options[0][1])))
        selection_numbers = [0]
    


    if args.output_filename != None:
        while len(selection_numbers) == 0: #continue prompt until a valid input is received
            user_in = input("Please select one or more sets of individuals to output (separate with spaces): ")
            inputs = user_in.split()
            for input_num in inputs:
                if not input_num.isdigit() or int() < 0 or int(input_num) >= len(list_options):
                    print("Invalid input. Please input one or more numbers from the lists above.")
                    selection_numbers = []
                    break
                else:
                    selection_numbers.append(int(input_num))

        for num in selection_numbers:
            output_list += list_options[num][1]
        #write chosen pedigree to output file
        write_to_file(args.output_filename,ped,list(set(output_list)),args.quiet)

    #create ped file
    if args.pedigree_filenames != None:
        create_ped_file(args.pedigree_filenames[0], args.pedigree_filenames[1], output_list,args.quiet)

def write_to_file(filename,ped,output_list,quiet):
    """
    writes ped struct to output file
    """
    out_file = open(filename, "w")
    out_file.write("ID FATHER MOTHER SEX")
    for id in output_list:
        indv = ped.indvs[id]
        #changes parents to 0 if they are not in the pedigree
        p = '0'
        m = '0'
        if indv.p_id in output_list and indv.m_id in output_list:
            p = indv.p_id
            m = indv.m_id
        out_file.write("\n" + id + " " \
            + p + " " \
            + m + " "  \
            + str(indv.sex))
    out_file.close()
    if not quiet:
        print("pedigree structure stored in " + filename)

def create_ped_file(input,output,ids,quiet):
    """
    copies only the minimum pedigree members
    from the input .ped file to the output .ped
    file
    """
    in_file = open(input, "r")
    lines = []

    #find relevant individuals in input file
    for line in in_file:
        if line.strip():
            words = line.split()
            if words[1] in ids:
                lines.append(line)
    in_file.close()

    #write relevant individuals to output file
    out_file = open(output, "w")
    for line in lines:
        out_file.write(line)
    out_file.close()

    if not quiet:
        print("pedigree contents stored in " + output)


def find_min_pedigree(ped,start_ids,source,timeout,quiet):
    """
    Takes a pedigree (pedigreeTree) and a list of ids (strings).
    Will find all shared ancestors for the starting indvs and get a list
    of minimum members to cover all paths from the source to the indvs.
    Returns a list of all found individual lists.
    """
    #print("starting ids: " + str(start_ids))

    #find shared ancestors
    shared_ancestors = ped.find_collective_ca(start_ids)
    #print("shared ancestors: " + str(len(shared_ancestors)))


    ped_options = []

    if source != None:
        source = source.replace("+","&")

    #create a list for each ancestor
    anc_count = 1
    for ancestor_id in shared_ancestors:

        if source != None and ancestor_id != source:
            continue

        if not quiet:
            print("finding set for ancestor " + str(anc_count)+"/"+str(len(shared_ancestors)), end='\r')

        #find all paths from ancestor to descendents
        ancestor = shared_ancestors[ancestor_id]

        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(timeout)
        try: #timeout if path finding takes too long
            all_paths = set() #set of ancestorNodes
            #get all paths from each start ids to source
            for id in start_ids:
                path_set = {}
                ped.get_all_paths(ancestor,id,ancestor.indv.sex,path_set)
                all_paths = all_paths | path_set.keys()
        except TimeoutException: #skip source at timeout
            if not quiet:
                print("timeout reached")
            anc_count += 1
            continue
        signal.alarm(0)

        min_ids = []
        for node in all_paths:
            ids = node.indv.id.split('&')
            for id in ids:
                min_ids.append(id)
            indv = ped.indvs[id]
            #add parents
            if node.indv.id != ancestor_id and indv.p != None and indv.m != None:
                min_ids.append(indv.m_id)
                min_ids.append(indv.p_id)
        min_ids = list(set(min_ids))

        ped_options.append((ancestor_id,min_ids))

        anc_count += 1
    return ped_options

    

if __name__ == "__main__":
    main()
