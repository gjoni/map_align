# _map_align_

_map_align_ takes two contact maps and returns an alignment that attempts to maximize the number of overlapping contacts while minimizing the number of gaps [1].

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

##### More details on some options:

* `-M MAX` templates with more than `MAX` residues will be skipped 
 (longer templates could much more time to be aligned)

* `-T TM` top N partial hits will be cleaned to exclude structurally similar matches: 
 if two hits from the scan stage are similar with TM-score > `TM` then only one 
 (with the higher alignment score) will appear in the final pool

* `-O PREFIX` top models will be saved as `<PREFIX><TEMPLATE_ID>.pdb` where `TEMPLATE_ID` is an ID from the `list` file


### Examples

Align a contact map to a PDB file (with and without saving the partial thread):
```
$ ./map_align -s example/T0806.fas -c example/T0806.con -p example/1A1X_A.pdb
$ ./map_align -s example/T0806.fas -c example/T0806.con -p example/1A1X_A.pdb -o example/T0806.1A1X_A.pdb
```

Align a contact map to a library of templates saving top 5 hits at TM-score=70% identity cut-off and running the program on 4 cores:
```
$ ./map_align -s example/T0806.fas -c example/T0806.con -D example -L example/list -N 5 -T 0.70 -t 4 -O example/T0806.
```


### Acknowledgements

This package is a reimplementation of the orinal map_align program by S.Ovchinnikov [1] https://github.com/sokrypton/map_align to allow for:
 - direct use of PDB files as templates
 - output of partial threads in PDB format
 - cleaning of partial threads based on TM-score [2]
 - multithreading

External packages/libraries:
 - kdtree library by John Tsiombikas https://github.com/jtsiomb/kdtree
 - C++ TM-align routine from Yang Zhang lab https://zhanglab.ccmb.med.umich.edu/TM-align

### References

[1] S Ovchinnikov et al. Protein Structure Determination using Metagenome sequence data. (2017) Science. 355(6322):294–8.

[2] Y Zhang & J Skolnick. TM-align: a protein structure alignment algorithm based on the TM-score. (2015) Nucleic Acids Res. 33(7):2302-9.

