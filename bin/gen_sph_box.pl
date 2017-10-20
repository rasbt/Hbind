#!perl
#
# -- Rajesh Korde Wed Apr 17 11:45:03 EDT 2002
# Usage: gen_sph_box.pl <x> <y> <z> <r>
# make a cube around a sphere
#

# Usage check
if (@ARGV != 4) {
  die "Usage: $0 <x> <y> <z> <r>\n";
}

$x=$ARGV[0];
$y=$ARGV[1];
$z=$ARGV[2];
$r=$ARGV[3];

printf("%8.2f %8.2f %8.2f\n", $x-$r, $y-$r, $z-$r);
printf("%8.2f %8.2f %8.2f\n", $x+$r, $y-$r, $z-$r);
printf("%8.2f %8.2f %8.2f\n", $x-$r, $y+$r, $z-$r);
printf("%8.2f %8.2f %8.2f\n", $x+$r, $y+$r, $z-$r);
printf("%8.2f %8.2f %8.2f\n", $x-$r, $y-$r, $z+$r);
printf("%8.2f %8.2f %8.2f\n", $x+$r, $y-$r, $z+$r);
printf("%8.2f %8.2f %8.2f\n", $x-$r, $y+$r, $z+$r);
printf("%8.2f %8.2f %8.2f\n", $x+$r, $y+$r, $z+$r);
