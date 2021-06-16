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
import pickle
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

class SubPedigree:
    def __init__(self,source,cohorts,mem_ids):
        self.source = source
        self.cohorts = cohorts
        self.mem_ids = mem_ids

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
    parser.add_argument("-pikl", "--pickle_filename", \
        help="a pickle file for saving found subpeds")
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

    #test_code(ped,IBDs,args,timeout)
    alternative_approach(ped,IBDs,args,timeout)

    
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
            print("[" + str(i) + "] source: " + list_options[i][0] + " pedigree size: " + str(len(list_options[i][1])) + " has loops: " + str(list_options[i][2]))
    else: #using preselected source
        print("source: " + list_options[0][0] + " pedigree size: " + str(len(list_options[0][1])) + " has loops: " + str(list_options[0][2]))
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

        contains_loops = False

        min_ids = []
        parent_ids = []
        for node in all_paths:
            ids = node.indv.id.split('&')
            for id in ids:
                min_ids.append(id)
            indv = ped.indvs[id]
            #add parents
            if node.indv.id != ancestor_id and indv.p != None and indv.m != None:
                parent_ids.append(indv.m_id)
                parent_ids.append(indv.p_id)
        parent_ids = list(set(parent_ids))

        married_in_count = 0
        for id in parent_ids:
            if not id in min_ids:
                married_in_count += 1
        
        for id in min_ids:
            indv = ped.indvs[id]
            if indv.p != None and indv.m != None \
                and indv.p_id + "&" + indv.m_id != ancestor_id \
                and indv.p_id in min_ids and indv.m_id in min_ids:
                contains_loops = True
                break

        min_ids = list(set(min_ids + parent_ids))

        subped = SubPedigree(ancestor_id,[start_ids],min_ids)

        ped_options.append(subped)

        anc_count += 1
    
    if not quiet:
        #print("                                                          ", end='\r')
        print("\033[K",end='\r')
    return ped_options

def alternative_approach(ped,IBDs,args,timeout):
    


    list_options = []
    if args.pickle_filename != None and os.path.exists(args.pickle_filename):
        pickle_file = open(args.pickle_filename,"rb")
        list_options = pickle.load(pickle_file)
    else:
        for selected_ibd in IBDs:
            starting_indvs = []
            for indv in selected_ibd.get_indvs():
                starting_indvs.append(indv)
            list_options += find_min_pedigree(ped,starting_indvs,args.source,timeout,args.quiet)
        if args.pickle_filename != None:
            pickle_file = open(args.pickle_filename,"wb")
            pickle.dump(list_options,pickle_file)
            pickle_file.close()

    source_options = {}
    for option in list_options:
        numeric_id = float(option.source.replace('&','.'))
        if numeric_id in source_options:
            source_options[numeric_id] += [option]
        else:
            source_options[numeric_id] = [option]


    #TODO change this to remove equivalent peds
    print("removing redundant peds...",end='\r')
    for source in source_options.keys():
        options = source_options[source]
        new_options = []
        for option in options:
            if not option in new_options:
                new_options.append(option)
        source_options[source] = new_options
    print("\033[K",end='\r')
    '''
    equal_peds = 0
    for source in source_options.keys():
        options = source_options[source]
        new_options = []
        for i in range(len(options)):
            for j in range(i+1,len(options)):
                if set(options[i].mem_ids) == set(options[j].mem_ids):
                    equal_peds += 1
    print(equal_peds)
    '''

    i = 0
    sorted_ids = []
    for id in sorted(source_options.keys()):
        if i % 20 == 0:
            print("source\t\t\tped count\ttotal mems")
        options = source_options[id]
        full_ped = []
        for option in options:
            full_ped += option.mem_ids
        full_ped = list(set(full_ped))

        source_id = options[0].source

        sorted_ids.append(source_id)
        spacing = ""
        while len("[" + str(i) + "] " + source_id + spacing) < 17:
            spacing += " "
        print("[" + str(i) + "] " + source_id + spacing + "\t" + str(len(options)) + "\t\t" + str(len(full_ped)))
        i+= 1
    
    selected_source = None
    while selected_source == None: #continue prompt until a valid input is received
        user_in = input("Please select a source from the list above: ")
        # TODO make multiple selections?
        if not user_in.isdigit() or int(user_in) < 0 or int(user_in) >= len(sorted_ids):
            print("Invalid input. Please input a number from the lists above.")
        else:
            selected_source = sorted_ids[int(user_in)]

    numeric_id = float(selected_source.replace('&','.'))
    options = source_options[numeric_id]

    full_ped = []
    subpeds =[]
    min_size = len(options[0].mem_ids)
    for option in options:
        full_ped += option.mem_ids
        subpeds.append(option)
        if len(option.mem_ids) < min_size:
            min_size = len(option.mem_ids)
    full_ped = list(set(full_ped))
    print("combined ped sizes range from " + str(min_size) + " to " + str(len(full_ped)))

    joined_ped = None
    while joined_ped == None:
        user_in = input("Please select a desired ped size: ")
        if not user_in.isdigit() or int(user_in) < min_size or int(user_in) >= len(full_ped):
            print("Invalid input. Please input a number in the range above.")
        else:
            target_size = int(user_in)

            table = []
            for i in range(len(subpeds)):
                row = []
                for j in range(len(full_ped)):
                    row.append(None)
                table.append(row)

            if sys.getrecursionlimit() < len(full_ped)**2:
                sys.setrecursionlimit(len(full_ped)**2) 

            low_option,high_option = join_peds(table,subpeds,target_size,SubPedigree(selected_source,[],[]),0)
            if low_option == high_option:
                joined_ped = low_option
                print("found ped of exact size")
                #print(table)
            else:
                print("could not be matched exactly, closest sizes are " + str(len(low_option.mem_ids)) + " and " + str(len(high_option.mem_ids)))
    

    #create output file
    if args.output_filename != None:
        write_to_file(args.output_filename,ped,list(set(joined_ped.mem_ids)),args.quiet)

    #create ped file
    if args.pedigree_filenames != None:
        create_ped_file(args.pedigree_filenames[0], args.pedigree_filenames[1], list(set(joined_ped.mem_ids)),args.quiet)




    exit()

def join_peds(table,subpeds,target_size,current_ped,list_num):
    '''
    returns closest subpeds below and above target size
    '''
    list_size = len(current_ped.mem_ids)

    if table[list_num][list_size] != None:
        return table[list_num][list_size][0],table[list_num][list_size][1]

    new_cohorts = current_ped.cohorts + subpeds[list_num].cohorts
    new_mem_ids = list(set(current_ped.mem_ids + subpeds[list_num].mem_ids))
    new_ped = SubPedigree(current_ped.source,new_cohorts,new_mem_ids)

    if len(new_mem_ids) == target_size:
        table[list_num][list_size] = (new_ped,new_ped)
        return new_ped,new_ped
    elif list_num + 1 == len(subpeds):
        table[list_num][list_size] = (current_ped,new_ped)
        return current_ped,new_ped
    
    low_option = current_ped
    high_option = None

    if len(new_ped.mem_ids) > target_size:
        high_option = new_ped

    else:
        low_option = new_ped
        recurse_low,recurse_high = join_peds(table,subpeds,target_size,new_ped,list_num+1)
        high_option = recurse_high
        if len(recurse_low.mem_ids) > len(low_option.mem_ids):
            print("rl2 chose " + str(len(recurse_low.mem_ids)) + " over " + str(len(low_option.mem_ids)))
            low_option = recurse_low

    if len(low_option.mem_ids) != target_size and len(high_option.mem_ids) != target_size:
        recurse_low,recurse_high = join_peds(table,subpeds,target_size,current_ped,list_num+1)
        if len(recurse_low.mem_ids) > len(low_option.mem_ids):
            print("rl1 chose " + str(len(recurse_low.mem_ids)) + " over " + str(len(low_option.mem_ids)))
            low_option = recurse_low
        if len(recurse_high.mem_ids) >= target_size and len(recurse_high.mem_ids) < len(high_option.mem_ids):
            print("rh1 chose " + str(len(recurse_high.mem_ids)) + " over " + str(len(high_option.mem_ids)))
            high_option = recurse_high
    
    if len(low_option.mem_ids) == target_size:
        high_option = low_option
    elif len(high_option.mem_ids) == target_size:
        low_option = high_option
    
    table[list_num][list_size] = (low_option,high_option)
    return low_option,high_option


    
    




def test_code(ped,IBDs,args,timeout):
    #testing code
    # python3 ped-cohort.py input/extended_pedigree_final.fam output/trial1_chr21_amr_geno_germline.match -o output/trial1_chr21_cohort.fam -p output/trial1_chr21_amr_geno.ped output/trial1_chr21_amr_cohort.ped -t 10
    
    
    i = 0
    #info_file = open("recon_info.txt","w")
    #info_file.close()
    #info_file = open("geno_info.txt","w")
    #info_file.close()
    #info_file = open("loop_info.txt","w")
    #info_file.close()

    for selected_ibd in IBDs:

        starting_indvs = []
        for indv in selected_ibd.get_indvs():
            starting_indvs.append(indv)
        list_options = find_min_pedigree(ped,starting_indvs,args.source,timeout,args.quiet)
        
        if False and args.output_filename != None:
            for option in list_options:
                write_to_file(args.output_filename,ped,option[1],args.quiet)
                status = os.system("Rscript visualizePed.R " + str(args.output_filename))
                if status != 0:
                    print("failed early with status: " + str(status))
                    exit()

        if False:
            info_file = open("loop_info.txt","a")
            info_file.write(selected_ibd.id + "\n")
            info_file.close()
            for option in list_options:
                info_file = open("loop_info.txt","a")
                info_file.write(str(option[0]) + " " + str(option[3]) + "\n")
                info_file.close()
        
        if False and args.pedigree_filenames != None:
            info_file = open("geno_info.txt","a")
            info_file.write(selected_ibd.id + "\n")
            info_file.close()
            for option in list_options:
                write_to_file(args.output_filename,ped,option[1],args.quiet)
                create_ped_file(args.pedigree_filenames[0], args.pedigree_filenames[1], option[1],args.quiet)
                ped_file = open(args.pedigree_filenames[1],"r")
                num_genotyped = 0
                for line in ped_file:
                    genotyped = False
                    if line.strip():
                        words = line.split()
                        for j in range(6,len(words)):
                            if words[j] != '?':
                                genotyped = True
                                break
                        if genotyped:
                            num_genotyped += 1
                info_file = open("geno_info.txt","a")
                info_file.write(str(num_genotyped) + "\n")
                info_file.close()
                

        if False and args.pedigree_filenames != None:
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
                        for j in range(6,len(words)):
                            if words[j] != '?':
                                didReconstruct = True
                                break
                        if didReconstruct:
                            break
                ped_file.close()
                info_file = open("recon_info.txt","a")
                info_file.write(option[0] + " "  + str(len(starting_indvs)) + " " + str(len(option[1])) + " " + str(option[2]) + " " + str(didReconstruct) +"\n")
                info_file.close()
        print("completed IBD " + str(i))
        i+=1
    
    print("all IBDs checked successfully")
    exit()
    

    

if __name__ == "__main__":
    main()
