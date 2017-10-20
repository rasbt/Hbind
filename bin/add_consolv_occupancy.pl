#!perl
#
# Time-stamp: <Paul Sanschagrin -- Thu Dec 20 11:00:37 EST 2001>
# add_consolv_occupancy.pl -- add Consolv results to a pdb file as 
#      the occupancy values for its waters as votes[conserved]/total 
#      number of votes Waters not predicted on will be set to 0.00.
# Usage: add_consolv_occupancy.pl [-a] <PDB Code>
#             -a: include all waters, default is only those with
#                 predictions > 0.50 (50%)
#
#

use File::Basename;
use Getopt::Std;

# The following is a hash to hold the valid water residue names
# These correspond to the water names in Consolv
%watnames = ("HOH" => 1,
	    "H2O" => 1,
	    "DOD" => 1,
	    "D20" => 1,
	    "WAT" => 1);
sub bynum {$a <=> $b}
sub usage () {
  $this = basename ($0);
  print "Usage: $this [-a] <PDB code>\n";
  print "       -a: include all waters, default is only those with\n";
  print "           predictions > 0.50 (50%)\n";
  print "*** This is run AFTER Consolv ***\n";
  exit 1;
}

# Usage check and CL arg parsing
usage() if (getopts('a') == false);
usage() unless @ARGV == 1;

$pdbcode = $ARGV[0];

# Set some filenames based on the PDB code
$pdbfile = "$pdbcode.pdb";
$hitsfile = "$pdbfile.hits";
$envfile = "$pdbcode.env";
$predfile = "$pdbcode.pred";
$outfile = "$pdbcode.predwats.pdb";

# check that we can open each of the necessary files
open PDB, $pdbfile or die "Unable to open PDB file: $pdbfile\n";
open HITS, $hitsfile or die "Unable to open hits file: $hitsfile\n";
open ENV, $envfile or die "Unable to open env file: $envfile\n";
open PRED, $predfile or die "Unable to open pred file: $predfile\n";
open OUT, ">$outfile" or die "Unable to open output file: $outfile\n";

# parse the env file (we need to do this since Consolv only predicts on
#   waters in the first shell, i.e. those with ADN > 0.0
$count = 0;
$numwats = 0;
while (<ENV>) {
  next if $_ =~ /^\#/; # skip comment lines
  s/^\s+//g,$_;
  ($code, $resnum, $adn) = split(/\s+/,$_);
  $adnarr{$resnum} = $adn;
  $newocc{$resnum} = 0.0;
  if ($adn > 0.0) { # here we can extract the index of the predicted water
    $resarr{$count} = $resnum;
    $count++;
  }
  $numwats ++;
}
close ENV;


# Now we must parse the prediction file
while (<PRED>) {
  next if $_ !~ /predicted/;
  ($junk,$index,$junk,$junk,$pred,$junk,$votesdisp, $votescons) = split(/\s+/,$_);
  if (exists ($resarr{$index})) {
    $newocc{$resarr{$index}} = $votescons/($votescons+$votesdisp);
  }
}

close PRED;

# Now we must parse the actual PDB file
while (<PDB>) {
  if ($_ !~ /^HETATM/) { # Can't be a water as it's not a HETATM
    print OUT $_;
    next;
  }
  $resname = substr($_,17,3);
  $resnum = substr($_,22,4)+0;
  if (! exists $watnames{$resname}) { # check if this HETATM is a water
    print OUT $_; # not a water, so print as is
    next;
  }
  # Now we've hit a waters and must output with the corrected occupancy
  if ($opt_a) { # we're outputting them all
    printf(OUT "%54s%6.2f%14s\n",substr($_,0,54),$newocc{$resnum},substr($_,60,14));
  } elsif ($newocc{$resnum} > 0.50) { # output only those predicted to be conserved
    printf(OUT "%54s%6.2f%14s\n",substr($_,0,54),$newocc{$resnum},substr($_,60,14));
  }
}

close PDB;
close OUT;
