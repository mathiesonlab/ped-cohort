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

Other options:

`-i [chromosome] [start] [end] [id.haplotype]` - Allows inputing a preselected IBD at the command line, bypassing the choice at runtime.

`-s [source]` - Allows inputing a preselected source at the command line, bypassing the choice at runtime. Sources should either formated as `id` an individual or `id1+id2` for a couple.

`-t [timeout]` - An integer number of seconds allowed for finding all paths from a given source to the selected IBD cohort. Defaults to 10 seconds.

`-q` - suppresses all terminal output except that needed for user input.

created by:
* Alton Wiggers (`ahwiggers`)

requires `thread` package.
