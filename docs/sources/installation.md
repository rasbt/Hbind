## Installation

Installing Hbind requires [Perl](https://www.perl.org) and the [GCC](https://gcc.gnu.org) compiler, both of which come pre-installed with most Unix- and Linux-based operating systems.

To install Hbind, simply download this repository, unpack it, and navigate into the main Hbind folder. Then, execute the following command in your terminal:

    perl install_hbind.pl

Upon successful installation, the Hbind software will be ready to use from the `bin/` subdirectory. To show the help and usage menu execute execute the following command in your terminal:

    ./bin/hbind -h

<br>

```
HBIND Version: 1.0.0

Documentation: http://psa-lab.github.io/Hbind
Raschka, Wolf, Bemister-Buffington, Kuhn (2018)
Protein Structure and Analysis Lab, MSU (http://kuhnlab.bmb.msu.edu)


USAGE:
-p STRING     Path to protein PDB file
-l STRING     Path to ligand mol2 file (in docked conformation)
-s            Include saltbridges in the output
-t            Include a summary table in the output
```


 Please see the "[User Guide](user_guide.md)" section for more information on how to use Hbind.