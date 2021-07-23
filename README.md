# ped-cohort

`ped-cohort` is a tool for finding minimal pedigrees that has only necessary path members from a source to all members of an IBD cohort.

genotype or sequence data should be phased (recommended: `SHAPEIT2`) and then IBD segments
should be called (recommended: `GERMLINE`). Then `ped-cohort` can be applied in the following way:

`ped-cohort` can be executed directly from the command line by running:

`$ python3 ped-cohort.py [struct file] [germline file] -o [output file]`

`struct file` should be a .txt file containing the pedigree structure, with each line being of the format:

`id father mother sex`

where `father` and `mother` are other ids in the file and sex is `1` for male, `2` for female, and `0` for unknown.

`germline file` should be a .match file produced by `GERMLINE`.

`output file` will be a text file in the same format as `struct file`.

`ped-cohort` will prompt the user to select a single IBD source identified in the pedigree. The user will then select a desired pedigree size and `ped-cohort` will merge minimum pedigrees containing various IBD cohorts and the selected source to get a pedigree of the desired size.

Other options:

`-p [input pedigree] [output pedigree]` - Allows for copying only the members in the chosen sub-pedigree from an associated `.ped` file to a new `.ped` file.
The `.ped` file should follow the plink ped format (https://www.cog-genomics.org/plink2/formats#ped).

`-c [component file name]` - A prefix for .ped files and corresponding .txt struct files for all component IBD cohort pedigrees that made up the final chosen pedigree. This will output `[component_filename]_[i].ped` and `[component_filename]_[i].txt` for the `i`th component. It's recommended to have these files output in a subdirectory because this can result in many files.

`-s [source]` - Allows inputing a preselected source at the command line, bypassing the choice at runtime. Sources should either formated as `id` for an individual or `id1+id2` for a couple.

`-m [maximum component complexity]` - Sets an integer number for a maximum bit complexity for component IBD cohort pedigrees. Bit complexity is `2n-f-g` where `n` is the number of non-founding pedigree members, `f` is the number of founding members, and `g` is the number of ungenotyped founding couples.

`pikl [pickle file]` - Allows the use of the python `pickle` package to save all sources and their assigned IBDs. If the file name does not exist in the current directory, `ped-cohort` will save a pickle file under this name. If the file does exist, `ped-cohort` will load the given pickle file.

`-q` - suppresses all terminal output except that needed for user input.

created by:
* Alton Wiggers (`ahwiggers`)

uses `AncestorNode.py`, `Couple.py`, `IBD.py`, `Individual.py`, and `PedigreeTree.py` from `thread` package (https://github.com/mathiesonlab/thread).
