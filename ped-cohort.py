"""
Main file to find minimum pedigree from an IBD cohort to a source
data.
Authors: Alton Wiggers
Date: 7/8/21
"""

#python imports
import argparse
import sys
import pickle
import os #used for testing

#local imports
import IBD
from PedigreeTree import PedigreeTree
from AncestorNode import AncestorNode

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
        help="name for output file.")
    parser.add_argument("-p", "--pedigree_filenames", nargs=2, \
        help="input and output file names for an associated .ped file (2 arg format: -p [input_ped] [output_ped]")
    parser.add_argument("-c", "--component_filename", \
        help="a prefix for .ped files for all components that made up the final pedigree will output [component_filename]_i.ped for i components.")
    parser.add_argument("-s", "--source", \
        help="a specific source to choose. Please use + instead of & for couples.")
    parser.add_argument("-m", "--max_component_size", type=int, \
        help="the maximum bit complexity for sub-pedigrees to consider when joining sub-pedigrees to reach a target size")
    parser.add_argument("-pikl", "--pickle_filename", \
        help="a pickle file for saving found subpeds")
    parser.add_argument("-q", "--quiet", action="store_true", \
        help="supress terminal output")

    args = parser.parse_args()

    return args


def main():

    #parse arguments
    args = parse_args("pedigree args")

    # construct pedigree data structure
    ped_tree = PedigreeTree(args.struct_filename)
    left_out = []

    # construct IBDs data structures (type: list)
    IBDs = IBD.get_IBDs(args.germ_filename, left_out) # toggle to leave out

    # assign IBDs to individuals in pedigree
    IBD.ibd_to_indvs(IBDs, ped_tree)

    #get a dictionary of source IDs and their minimum possible subpedigrees
    source_options = get_source_options(ped_tree,IBDs,args)

    #let user select a source and desired pedigree size
    chosen_ped = get_user_selection(ped_tree,args,source_options)

    #create output file
    if args.output_filename != None:
        write_to_file(args.output_filename,ped_tree,list(set(chosen_ped.mem_ids)),args.quiet)

    #create ped file
    if args.pedigree_filenames != None:
        create_ped_file(args.pedigree_filenames[0], args.pedigree_filenames[1], list(set(chosen_ped.mem_ids)),args.quiet)


def write_to_file(filename,ped_tree,output_list,quiet):
    """
    writes ped struct to output file
    """
    out_file = open(filename, "w")
    out_file.write("ID FATHER MOTHER SEX")
    for id in output_list:
        indv = ped_tree.indvs[id]
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


def find_min_pedigree(ped_tree,start_ids,source,quiet):
    """
    Takes a pedigree (pedigreeTree) and a list of ids (strings).
    Will find all shared sources for the starting indvs and get a SubPedigree
    of minimum members to cover all paths from the source to the indvs.
    Returns a list of all found SubPedigrees.
    """
    #print("starting ids: " + str(start_ids))

    #find shared ancestors
    shared_ancestors = ped_tree.find_collective_ca(start_ids)

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

        all_paths = set() #set of ancestorNodes
        #get all paths from each start ids to source
        for id in start_ids:
            path_set = {}
            ped_tree.get_all_paths(ancestor,id,ancestor.indv.sex,path_set)
            all_paths = all_paths | path_set.keys()

        contains_loops = False

        min_ids = []
        parent_ids = []
        for node in all_paths:
            ids = node.indv.id.split('&')
            for id in ids:
                min_ids.append(id)
            indv = ped_tree.indvs[id]
            #add parents
            if node.indv.id != ancestor_id and indv.p != None and indv.m != None:
                parent_ids.append(indv.m_id)
                parent_ids.append(indv.p_id)
        parent_ids = list(set(parent_ids))

        married_in_count = 0
        for id in parent_ids:
            if not id in min_ids:
                married_in_count += 1
        
        for id in min_ids: #identify loops in the pedigree
            indv = ped_tree.indvs[id]
            if indv.p != None and indv.m != None \
                and indv.p_id + "&" + indv.m_id != ancestor_id \
                and indv.p_id in min_ids and indv.m_id in min_ids:
                contains_loops = True
                break
        #sort ids for later comparisons
        min_ids = sorted(list(set(min_ids + parent_ids)))

        subped = SubPedigree(ancestor_id,[start_ids],min_ids)

        ped_options.append(subped)

        anc_count += 1
    
    if not quiet:
        print("\033[K",end='\r')
    return ped_options

def get_source_options(ped_tree,IBDs,args):
    """
    get a dictionary of sources and a list of SubPedigree
    objects for each IBD cohort with the keyed source as
    a possible source.
    """

    list_options = []
    #check for pickle file
    if args.pickle_filename != None and os.path.exists(args.pickle_filename):
        pickle_file = open(args.pickle_filename,"rb")
        list_options = pickle.load(pickle_file)
    else:
        #get options from each IBD cohort
        for selected_ibd in IBDs:
            starting_indvs = []
            for indv in selected_ibd.get_indvs():
                starting_indvs.append(indv)
            list_options += find_min_pedigree(ped_tree,starting_indvs,args.source,args.quiet)
        #save to pickle file
        if args.pickle_filename != None:
            pickle_file = open(args.pickle_filename,"wb")
            pickle.dump(list_options,pickle_file)
            pickle_file.close()

    source_options = {}
    #traverse all options and assign them to the correct source
    for option in list_options:
        if option.source in source_options:
            source_options[option.source] += [option]
        else:
            source_options[option.source] = [option]

    if not args.quiet:
        print("removing redundant peds...",end='\r')
    #remove redundant pedigrees
    for source in source_options.keys():
        options = source_options[source]
        new_options = []
        mem_lists = []
        #compares mem_lists of each SubPedigree
        for option in options:
            if not option.mem_ids in mem_lists:
                new_options.append(option)
                mem_lists.append(option.mem_ids)
        source_options[source] = new_options
    if not args.quiet:
        print("\033[K",end='\r')

    return source_options
    
def get_user_selection(ped_tree,args,source_options):
    """
    Prompts user for source selection and target pedigree
    size. Returns a SubPedigree of the proper size.
    """

    selected_source = None

    if args.source != None: #use preselected source if given in command line
        selected_source = args.source.replace('+','&')
    else:
        i = 0
        sorted_ids = []
        #print an option for each source
        for id in sorted(source_options.keys()):
            options = source_options[id]
            full_ped = []
            valid_options = 0
            #find the union of valid subpeds given maximum allowed complexity
            for option in options:
                if args.max_component_size == None or get_bit_complexity(ped_tree,option.mem_ids) <= args.max_component_size:
                    full_ped += option.mem_ids
                    valid_options += 1
            full_ped = list(set(full_ped))

            source_id = options[0].source

            if valid_options > 0: #only show options that have some valid subped
                if i % 20 == 0:
                    print("source\t\t\tped count\ttotal mems")
                sorted_ids.append(source_id)
                spacing = "" #adjusts spacing for clarity
                while len("[" + str(i) + "] " + source_id + spacing) < 17:
                    spacing += " "
                print("[" + str(i) + "] " + source_id + spacing + "\t" + str(valid_options) + "\t\t" + str(len(full_ped)))
                i+= 1
    
    while selected_source == None: #continue prompt until a valid input is received
        user_in = input("Please select a source from the list above: ")
        # TODO make multiple selections?

        #prevent invalid input
        if not user_in.isdigit() or int(user_in) < 0 or int(user_in) >= len(sorted_ids):
            print("Invalid input. Please input a number from the lists above.")
        else:
            selected_source = sorted_ids[int(user_in)]

    if not selected_source in source_options.keys():
        print("could not find selected source")
        exit()
    
    #get SubPedigrees for the chosen source
    options = source_options[selected_source]

    full_ped = []
    subpeds =[]
    min_size = len(options[0].mem_ids)
    #find minimum and maximum size options
    for option in options:
        #only use valid SubPedigrees
        if args.max_component_size == None or get_bit_complexity(ped_tree,option.mem_ids) <= args.max_component_size:
            full_ped += option.mem_ids
            subpeds.append(option)
            #use smallest SubPedigree for minimum size
            if len(option.mem_ids) < min_size:
                min_size = len(option.mem_ids)
    full_ped = list(set(full_ped)) #use union of all SubPedigrees for maximum size
    print("for " + selected_source + " combined pedigree sizes range from " + str(min_size) + " to " + str(len(full_ped)))

    joined_ped = None
    #find SubPedigree based on user selected size
    while joined_ped == None:
        user_in = input("Please select a desired pedigree size in the range above: ")
        #prevent invalid inputs
        if not user_in.isdigit() or int(user_in) < min_size or int(user_in) > len(full_ped):
            print("Invalid input. Please input a number in the range above.")
        else:
            target_size = int(user_in)
            #search for options of target size
            low_option,high_option = find_joined_ped(selected_source,subpeds,target_size,len(full_ped))
            
            #found SubPedigree of exact specified size
            if low_option == high_option:
                joined_ped = low_option

                cohort_min = target_size
                cohort_max = 0
                #find minimum and maximum sizes of joined subpeds
                for min_ped in subpeds:
                    if min_ped.cohorts[0] in joined_ped.cohorts:
                        if len(min_ped.mem_ids) < cohort_min:
                            cohort_min = len(min_ped.mem_ids)
                        if len(min_ped.mem_ids) > cohort_max:
                            cohort_max = len(min_ped.mem_ids)

                #print outcome
                if not args.quiet:
                    print("found pedigree of desired size by joining sub-peds of sizes "  + str(cohort_min) + "-" + str(cohort_max)   + \
                    " from " + str(len(joined_ped.cohorts)) + " different IBD cohorts")
                
                if args.component_filename != None:
                    create_component_files(ped_tree,args,joined_ped,subpeds)
                
            #If not exact pedigree was found, show closest sizes and reprompt
            else:
                print("could not be matched exactly, closest sizes are " + str(len(low_option.mem_ids)) + " and " + str(len(high_option.mem_ids)))
    
    return joined_ped


def find_joined_ped(source,subpeds,target_size,max_size):
    """
    Base call for recursive algorithm to join subpeds to
    get a subped of target size. Uses dynamic approach.

    Uses a dynamic approach with a table that has # subpeds
    rows and maximum pedigree size columns (size of the union
    of all subpeds).

    Returns closest Subpedigrees above and below the
    target size.
    """
    #create table for dynamic approach
    table = []
    for i in range(len(subpeds)): #create rows
        row = []
        for j in range(max_size): #create columns
            row.append(None) #all cell start as 'None'
        table.append(row)

    #increase recursion limit
    if sys.getrecursionlimit() < max_size**2:
        sys.setrecursionlimit(max_size**2)

    #make first recursive call.
    low_option,high_option = join_peds(table,subpeds,target_size,SubPedigree(source,[],[]),0)
    return low_option,high_option


def join_peds(table,subpeds,target_size,current_ped,list_num):
    """
    Recursive step for joining SubPedigrees.
    returns closest subpeds below and above target size.

    At each iterationtable indexed by list_num (the first
    list_num cohort subpeds that have been checked) and list_size
    (the size of the current given SubPedigree).
    """
    list_size = len(current_ped.mem_ids)

    '''Base Cases'''

    #case: if the best SubPedigree has already been computed
    if table[list_num][list_size] != None:
        return table[list_num][list_size][0],table[list_num][list_size][1]

    #add the list_num cohort ped to the current ped
    new_cohorts = current_ped.cohorts + subpeds[list_num].cohorts
    new_mem_ids = list(set(current_ped.mem_ids + subpeds[list_num].mem_ids))
    new_ped = SubPedigree(current_ped.source,new_cohorts,new_mem_ids)

    #case: we've found a subped of target size
    if len(new_mem_ids) == target_size:
        table[list_num][list_size] = (new_ped,new_ped)
        return new_ped,new_ped
    #we are at the last cohort subped
    elif list_num + 1 == len(subpeds):
        table[list_num][list_size] = (current_ped,new_ped)
        return current_ped,new_ped
    
    '''Recursive Cases'''

    low_option = current_ped
    high_option = None

    
    if len(new_ped.mem_ids) > target_size:
        high_option = new_ped
    else: #only need to join more peds if new_ped size is less than target size
        low_option = new_ped #target_size > new_ped size >= current_ped size
        #case: join new_ped with subpeds after list_num (keep subpeds[list_num] in the union)
        recurse_low,recurse_high = join_peds(table,subpeds,target_size,new_ped,list_num+1)
        high_option = recurse_high #take the only high option
        #pick the better low option
        if len(recurse_low.mem_ids) > len(low_option.mem_ids):
            low_option = recurse_low

    if len(low_option.mem_ids) != target_size and len(high_option.mem_ids) != target_size:
        #case: join current_ped with subpeds after list_num (don't keep subpeds[list_num] in the union)
        recurse_low,recurse_high = join_peds(table,subpeds,target_size,current_ped,list_num+1)
        #pick the better low option
        if len(recurse_low.mem_ids) > len(low_option.mem_ids):
            low_option = recurse_low
        #pick the better high option
        if len(recurse_high.mem_ids) >= target_size and len(recurse_high.mem_ids) < len(high_option.mem_ids):
            high_option = recurse_high
    
    #if either option is of target size, set both options to be the target size option
    if len(low_option.mem_ids) == target_size:
        high_option = low_option
    elif len(high_option.mem_ids) == target_size:
        low_option = high_option
    
    #fill out the table
    table[list_num][list_size] = (low_option,high_option)
    return low_option,high_option


def create_component_files(ped_tree,args,full_ped,subpeds):
    """
    Creates a .ped file for each component that made up
    the final chosen pedigree.
    """

    #get list of components
    components = get_ped_components(full_ped,subpeds)
    pedfile_lines = open(args.pedigree_filenames[0],"r").readlines()

    line_dict = {}

    #find the number of SNP markers in each line of the pedigree
    markers = len(pedfile_lines[0].split()) - 6

    #create a dictionary of all lines in the full pedigree
    for line in pedfile_lines:
        words = line.split()
        line_dict[words[1]] = words

    
    for i in range(len(components)): #iterate through the components
        if not args.quiet:
            print("creating component file " + str(i+1) + "/" + str(len(components)),end='\r')
        #create files
        component = components[i]
        outfile_name = args.component_filename + "_" + str(i) + ".ped"
        textfile_name = args.component_filename + "_" + str(i) + ".txt"
        outfile = open(outfile_name, "w")
        textfile = open(textfile_name, "w")
        textfile.write("ID FATHER MOTHER SEX") #write header


        for id in component.mem_ids: #write line for each member of the component pedigree
            indv = ped_tree.indvs[id]
            #only put parents for an individual if they are also in the component pedigree
            dad = "0"
            mom = "0"
            if indv.p_id in component.mem_ids and indv.m_id in component.mem_ids:
                dad = indv.p_id
                mom = indv.m_id
            #write line to .txt file
            out_line = id + " " + dad + " " + mom + " " + str(indv.sex)
            textfile.write("\n" + out_line)
            #write line to .ped file
            out_line = "1 " + out_line #+ " 0"
            if id in line_dict.keys(): #write haplotypes if known
                words = line_dict[id]
                for j in range(6,len(words)):
                    out_line += " " + words[j]
            else: #write zeros if haplotypes are unknown
                for j in range(markers):
                    out_line += " 0"
            outfile.write(out_line + "\n")
        outfile.close()
        textfile.close()
        if not args.quiet:
            print("\033[K",end='\r')
    
    return len(components)

    

def get_ped_components(full_ped,subpeds):
    """
    Searches through a list of subpeds and
    find which are components of the full ped
    based on the cohorts in the full ped.
    """
    components = []
    for cohort in full_ped.cohorts:
        for subped in subpeds:
            if subped.cohorts[0] == cohort:
                components.append(subped)
                break
    #print(len(components))
    return components



def get_bit_complexity(ped_tree,mem_ids):
    """
    calculate the bit complexity of a pedigree
    based on a list of member ids.
    """
    n = 0
    f = 0
    g = []
    for id in mem_ids:
        indv = ped_tree.indvs[id]
        #check if parents are in the pedigree
        if indv.p_id in mem_ids and indv.m_id in mem_ids:
            n += 1
        else:
            f +=1
            #add couples that are both founders
            for couple in indv.couples:
                if not couple in g \
                    and couple.p.id in mem_ids and couple.m.id in mem_ids \
                    and not couple.p.p_id and not couple.p.m_id in mem_ids \
                    and not couple.m.p_id and not couple.m.m_id in mem_ids:
                    g.append(couple)

    
    return 2*n-f-len(g)

if __name__ == "__main__":
    main()
