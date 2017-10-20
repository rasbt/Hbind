#!perl
#
# Time-stamp: <Paul Sanschagrin -- Wed Jul 11 15:43:34 EDT 2001>
# replace_score.pl: replace the score in a hbind mol2 output file
# Usage: replace_score.pl <file> <new score>
#

# Usage Check and CL parsing
die "Usage: $0 <file> <new score>\n" unless @ARGV == 2;
$fn = $ARGV[0];
$newscore = $ARGV[1];
die "ERROR: score should consist only of [+-.0123456789]\n" 
  if $newscore =~ /[^\+\.\-\d]/;

# Open file
open FILE, $fn or die "Unable to open file: $fn\n";

while (<FILE>) {
    if ($_ =~ Affinity score) 
    {
	$_ = "# Affinity score (kcal/mol): $newscore\n";
    }
    print $_;
    if ($_ =~ affiscore) 
    {
	$_ = "# affiscore (kcal/mol) : $newscore\n";
    }
    print $_;
}

close FILE;

