# _map_align_

_map_align_ takes two contact maps and returns an alignment that attempts to maximize the number of overlapping contacts while minimizing the number of gaps.

![example image](https://raw.githubusercontent.com/sokrypton/map_align/master/map_align_fig.png)

### Download and Installation
```sh
$ git clone https://github.com/gjoni/map_align
$ cd map_align
$ make
```

### Usage

```
Usage:   ./map_align [-option] [argument]

Options:  -s sequence.fas                - input, required
          -c contacts.txt                - input, required

          ***************** single template ****************
          -p template.pdb                - input, required
          -o match.pdb                   - output, optional

                                  OR                        
          ************** library of templates **************
          -D path to templates           - input, required
          -L list.txt with template IDs  - input, required
          -O prefix for saving top hits  - output, optional
          -N number of top hits to save    (10)
          -T TM-score cleaning cut-off     (0.80)
          -M max template size             (1000)

          ********************** misc **********************
          -t number of threads             (1)

```

### Examples

Align a contact map to a PDB file (with and without saving the partial thread):

```
$ ./map_align -s sequence.fas -c contacts.txt -p template.pdb
$ ./map_align -s sequence.fas -c contacts.txt -p template.pdb -o match.pdb
```

Align a contact map to a library of templates saving top 30 hits at TM-score=70% identity cut-off and running the program on 4 threads:
```
$ ./map_align -s sequence.fas -c contacts.txt -D PATH -L list -N 30 -T 0.70 -t 4 -O model_
```


### Experimental features

### Acknowledgements

This package is a reimplementation of the orinal map_align program by S.Ovchinnikov [1] to allow for:
 - direct use of PDB files as templates
 - output of partial threads in PDB format
 - cleaning of partial threads based on TM-score [2]
 - multithreading

External packages/libraries:
 - kdtree by John Tsiombikas https://github.com/jtsiomb/kdtree
 - C++ TM-align routine from Y.Zhang lab https://zhanglab.ccmb.med.umich.edu/TM-align

References:
[1] S Ovchinnikov et al. Protein Structure Determination using Metagenome sequence data. (2017) Science. 355(6322):294–8.
[2] Y Zhang & J Skolnick. TM-align: a protein structure alignment algorithm based on the TM-score. (2015) Nucleic Acids Res. 33(7):2302-9.



