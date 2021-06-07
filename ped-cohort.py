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

default_timeout = 10

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
        help="input and output file names for an associated .ped file (2 arg format: input_ped output_ped")
    parser.add_argument("-i", "--ibd", nargs=4, \
        help="a specific ibd to choose (4 arg format: chr start end indv)")
    parser.add_argument("-s", "--source", \
        help="a specific source to choose (recommended to use with -i). Please use + instead of & for couples.")
    parser.add_argument("-t", "--timeout", type=int, \
        help="number of seconds allowed to search for paths for an individual ancestor (default is 10)")
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

    '''
    i = 0
    for selected_ibd in IBDs:
        starting_indvs = []
        for indv in selected_ibd.get_indvs():
            starting_indvs.append(indv)
        list_options = find_min_pedigree(ped,starting_indvs,args.source,timeout,args.quiet)
        
        if args.output_filename != None:
            for option in list_options:
                write_to_file(args.output_filename,ped,option[1],args.quiet)
                #status = os.system("R --quiet --vanilla < visualizePed.R " + str(args.output_filename))
                status = os.system("Rscript visualizePed.R " + str(args.output_filename))
                if status != 0:
                    print("failed early with status: " + str(status))
                    exit()
        print("completed IBD " + str(i))
        i+=1
    print("all IBDs checked successfully")
    exit()
    '''

    
    selected_ibd = None
    if args.ibd == None: #check for ibd selection at command line
        for i in range(len(IBDs)): # print options
            print("[" + str(i) + "] " + IBDs[i].id)
        
        while selected_ibd == None: #continue prompt until a valid input is received
            user_in = input("Please select an IBD from the list above: ")
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
    
    selection_number = -1

    if len(list_options) == 0:
        print("no pedigree options were found")
        exit()
    elif args.source == None: #check for source selection at command line
        for i in range(len(list_options)): #print options
            print("[" + str(i) + "] source: " + list_options[i][0] + " pedigree size: " + str(len(list_options[i][1])))
    else: #using preselected source
        print("source: " + list_options[0][0] + " pedigree size: " + str(len(list_options[0][1])))
        selection_number = 0
    


    if args.output_filename != None:
        while selection_number < 0: #continue prompt until a valid input is received
            user_in = input("Please select a set of individuals to output: ")
            if not user_in.isdigit() or int(user_in) < 0 or int(user_in) >= len(list_options):
                print("Invalid input. Please input a number from the lists above.")
            else:
                selection_number = int(user_in)

        output_list = list_options[selection_number][1]
        #write chosen pedigree to output file
        write_to_file(args.output_filename,ped,output_list,args.quiet)

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
    out_file = open(output, "w")

    for line in in_file:
        if line.strip():
            words = line.split()
            if words[1] in ids:
                out_file.write(line)
    
    in_file.close()
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
        #print(ancestor_id)
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(timeout)
        try: #timeout if path finding takes too long
            ancestor_paths = ped.descendence_paths({ancestor_id : ancestor }, start_ids)
        except TimeoutException:
            if not quiet:
                print("timeout reached")
            anc_count += 1
            continue
        signal.alarm(0)
        paths = ancestor_paths[ancestor]
        #print("found " + str(len(paths)) +  " paths")
        joined_paths = set().union(*paths)
        #print("all path mems: " + str(len(joined_paths)))


        #split ids for couples
        min_ids = ancestor.indv.id.split('&')
        #add members from each node in each path
        for node in joined_paths:
            indv = node[0].indv
            min_ids.append(indv.id)
            #add parents
            if indv.p != None and indv.m != None:
                min_ids.append(indv.m_id)
                min_ids.append(indv.p_id)
        
        #removes redundant members
        min_ids = list(set(min_ids))
        ped_options.append((ancestor_id,min_ids))
        #print(min_ids)
        anc_count += 1
    return ped_options

    

if __name__ == "__main__":
    main()
