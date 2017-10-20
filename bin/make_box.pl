#!perl
#
# Time-stamp: <Paul Sanschagrin -- Thu Aug  2 10:31:49 EDT 2001>
# make_box.pl: Make a box in a PDB formatted file given a set of 6 corners
# Usage: make_box.pl <x1> <y1> <z1> <x2> <y2> <z2>

# Usage check
if (@ARGV != 6) {
  die "Usage: $0 <x1> <y1> <z1> <x2> <y2> <z2>\n";
}

# Set some referencing keys
$x = 0;
$y = 1;
$z = 2;

# Output the box
printf ("HETATM    1  O   XXX     1    %8.3f%8.3f%8.3f\n",
	$ARGV[$x+0], $ARGV[$y+0], $ARGV[$z+0]);
printf ("HETATM    2  O   XXX     2    %8.3f%8.3f%8.3f\n",
	$ARGV[$x+3], $ARGV[$y+0], $ARGV[$z+0]);
printf ("HETATM    3  O   XXX     3    %8.3f%8.3f%8.3f\n",
	$ARGV[$x+0], $ARGV[$y+3], $ARGV[$z+0]);
printf ("HETATM    4  O   XXX     4    %8.3f%8.3f%8.3f\n",
	$ARGV[$x+3], $ARGV[$y+3], $ARGV[$z+0]);
printf ("HETATM    5  O   XXX     5    %8.3f%8.3f%8.3f\n",
	$ARGV[$x+0], $ARGV[$y+0], $ARGV[$z+3]);
printf ("HETATM    6  O   XXX     6    %8.3f%8.3f%8.3f\n",
	$ARGV[$x+3], $ARGV[$y+0], $ARGV[$z+3]);
printf ("HETATM    7  O   XXX     7    %8.3f%8.3f%8.3f\n",
	$ARGV[$x+0], $ARGV[$y+3], $ARGV[$z+3]);
printf ("HETATM    8  O   XXX     8    %8.3f%8.3f%8.3f\n",
	$ARGV[$x+3], $ARGV[$y+3], $ARGV[$z+3]);
print "CONECT    1    2\n";
print "CONECT    1    3\n";
print "CONECT    1    4\n";
print "CONECT    1    5\n";
print "CONECT    1    6\n";
print "CONECT    1    7\n";
print "CONECT    2    3\n";
print "CONECT    2    4\n";
print "CONECT    2    5\n";
print "CONECT    2    6\n";
print "CONECT    2    8\n";
print "CONECT    3    4\n";
print "CONECT    3    5\n";
print "CONECT    3    7\n";
print "CONECT    3    8\n";
print "CONECT    4    6\n";
print "CONECT    4    7\n";
print "CONECT    4    8\n";
print "CONECT    5    6\n";
print "CONECT    5    7\n";
print "CONECT    5    8\n";
print "CONECT    6    7\n";
print "CONECT    6    8\n";
print "CONECT    7    8\n";

