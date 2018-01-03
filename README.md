# Hbind -- Identifying hydrogen bonds by donor/acceptor chemistry and geometric constraints

### in progress ...

Software to rigorously define intermolecular H-bonds by donor/acceptor chemistry and geometric constraints, which was developed, used, and described in detail in 

- Sebastian Raschka, Wolf A., Bemister-Buffington J., and Kuhn L.A. (2018) 
"Protein-ligand interfaces are polarized: Discovery of a strong trend for intermolecular hydrogen bonds to favor donors on the protein side with implications for predicting and designing ligand complexes." J. Computer-Aided Molec. Design [in revision]

Documentation: [insert link]

## Installation

Installing Hbind requires [Perl](https://www.perl.org) and the [GCC](https://gcc.gnu.org) compiler, both of which come pre-installed with most Unix- and Linux-based operating systems.

To install Hbind, simply download this repository, unpack it, and navigate into the main Hbind folder. Then, execute the following command in your terminal:

    perl install_hbind.pl

Upon successful installation, the hbind software will be ready to use from the `bin/` subdirectory. To show the help menu execute

    ./bin/hbind_score --help

and to show a brief usage example, run

    ./bin/hbind_score --usage

in your command line terminal. Please see the "Usage" section below for more information on how to use Hbind.

## Usage

in progress ... 

## Hydrogen-bond rules

Hbind's rules for identifying hydrogen bonds protein-ligand interfaces are based on the criteria by Ippolito et al. [1] and McDonald and Thornton [2].

#### Hydrogen bond criteria:

- Hydrogen to acceptor distance: 1.5-2.5 Å
- Donor to acceptor distance: 2.4-3.5 Å
- Donor-H-acceptor angle (θ): 120-180°
- Pre-acceptor–acceptor–H angle (φ): 90-180°

(Hydrogen bonds must meet all of these criteria.)

#### Criteria for metal interactions:

- Lone pair atom distance to K or Na: 2.0-2.9 Å
- Lone pair atom distance to Ca, Co, Cu, Fe, Mg, Mn, Ni, or Zn: 1.7-2.6 Å

Additional command line options are available to list longer-range salt bridge interactions (up to 4.5 Å) and the number of hydrophobic contancts between protein and ligand atoms.


[1] Ippolito, Joseph A, Richard S Alexander & David W Christianson. 1990. Hydrogen bond stereo-chemistry in protein structure and function. *Journal of Molecular Biology* 215(3). 457–471.   
[2] McDonald, Ian & Janet M Thornton. 1994. Atlas of side-chain and main-chain hydrogen bond- ing. http://www.biochem.ucl.ac.uk/bsm/atlas: Biochemistry and Molecular Biology Department, University College London.

## Hbind tools

###  Protein Recognition Index (PRI) 

Software for assessing the native-likeness of designed or predicted protein-ligand interfaces, which can be used to guide protein mutagenesis and ligand design.

- [Link]

### Hbind Interaction Viz [insert link]

A program for creating PyMOL H-bond interaction views from Hbind tables.

- [Link]

