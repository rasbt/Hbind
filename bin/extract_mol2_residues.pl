#!perl
#
# extract_mol2_residues.pl  Volker Schnecke   Tue Feb 24 16:58:16 EST 1998
#
# Usage: extract_mol2_residues.pl directory
#
# This script reads all mol2 files in <directory> and checks if there
# are more than one residue in a file. If so, these are extracted 
# and stored in files <basename>_1.mol2, <basename>_2.mol2, ...

if ( $#ARGV != 0 )
{
    die "Usage: $0 <directory>\n";
}
 
opendir DIR, $ARGV[0];
@files = grep ( /.mol2$/, grep ( !/^.\.?$/, readdir ( DIR ) ) );
close DIR;

foreach $file ( @files )
{
    $file =~ s/.mol2//;

    open IN, "$ARGV[0]/$file.mol2";
    $_ = <IN>;	        # this should be "@<TRIPOS>MOLECULE"
    $compound = <IN>;	# read line containing compound name
    $_ = <IN>;		# this line includes numbers of atoms and bonds
    chop;
    @line = split;
    $number_of_atoms = $line[0];
    $number_of_bonds = $line[1];
    $number_of_residues = $line[2];

    if ( $number_of_residues == 1 )
    # nothing to do, so exit
    {
	close IN;
	next;
    }

    # read all lines up to the beginning of the atom definition and store
    # them in an array, this header will be included in each of the new files
    # note that the last line in this array is "@<TRIPOS>ATOM"
    $number_of_header_lines = 0;
    $header[$number_of_header_lines] = <IN>;
    while ( not $header[$number_of_header_lines] =~ /^@<TRIPOS>ATOM/ )
    {				
	$number_of_header_lines++;
	$header[$number_of_header_lines] = <IN>;
    }

    for ( $i = 1; $i <= $number_of_residues; $i++ )
    {
	$number_of_residue_atoms[$i] = 0;
	$number_of_residue_bonds[$i] = 0;
    }

    # read atom definitions
    # arrays from 1..n, since atoms and bonds are numbered this way
    for ( $i = 1; $i <= $number_of_atoms; $i++ )
    {
	$atoms[$i] = <IN>;
	$atoms[$i] =~ s/^\s+//;		# remove leading spaces
	@line = split /\s+/, $atoms[$i];
	$residue = $line[6];
	$atom_residue[$i] = $residue;
	# atom_index contains the new index of this atom in the 
	# particular residue
	$atom_index[$i] = $number_of_residue_atoms[$residue] + 1;
	$number_of_residue_atoms[$residue]++;
    }
    
    $_ = <IN>;			# this should be "@<TRIPOS>BOND", skip it
    # read bond definitions
    for ( $i = 1; $i <= $number_of_bonds; $i++ )
    {
	$bonds[$i] = <IN>;
	$bonds[$i] =~ s/^\s+//;
	@line = split /\s+/, $bonds[$i];
	$residue = $atom_residue[$line[1]];
	# new indices for both atoms adjacent to this bond
	$atom1[$i] = $atom_index[$line[1]];
	$atom2[$i] = $atom_index[$line[2]];
	$bond_type[$i] = $line[3];
	$bond_residue[$i] = $residue;
	$bond_index[$i] = $number_of_residue_bonds[$residue] + 1;
	$number_of_residue_bonds[$residue]++;
    }
    
    # skip the substructure part
    while ( <IN> )
    {
	last if /^@<TRIPOS>CRYSIN/;
    }

    # read the remaining lines
    $number_of_tail_lines = 1;
    $tail[0] = $_;
    while ( <IN> )
    {
	$tail[$number_of_tail_lines] = $_;
	$number_of_tail_lines++;
    }

    close IN;

    print "creating files ";
    for ( $res = 1; $res <= $number_of_residues; $res++ )
    {
      print "$file\_$res.mol2 ";
      # open file <basename>_$res.mol2 for this residue
      open OUT, ">$ARGV[0]/$file\_$res.mol2" or 
	die "Can't open output file: $ARGV[0]/$file\_$res.mol2\n";
      # write the header
      print OUT "@<TRIPOS>MOLECULE\n$compound";
      printf OUT "%6d %5d     1\n", 
	$number_of_residue_atoms[$res],
	  $number_of_residue_bonds[$res];
      for ( $i = 0; $i <= $number_of_header_lines; $i++ )
	{			
	  print OUT $header[$i];
	}
      # write the atom definitions
      for ( $i = 1; $i <= $number_of_atoms; $i++ )
	{
	  if ( $atom_residue[$i] == $res )
	    # only those atoms that are in this residue
	    {
	      $atoms[$i] =~ s/.\d RES\d./ 1 RES1 /;
	      $atoms[$i] =~ s/\s*\d*//;
	      printf OUT "%6d%s", $atom_index[$i], $atoms[$i];
	    }
	}
      # write the bond definitions
      print OUT "@<TRIPOS>BOND\n";
      for ( $i = 1; $i <= $number_of_bonds; $i++ )
	{
	  if ( $bond_residue[$i] == $res )
	    {
	      printf OUT "%6d%6d%6d    %s\n", 
		$bond_index[$i], 
		$atom1[$i],
		$atom2[$i],
		$bond_type[$i];
	    }
	}
      # of course, only one substructure
      print OUT "@<TRIPOS>SUBSTRUCTURE\n    1 RES1       1\n";
      # write tail of the file
      for ( $i = 0; $i < $number_of_tail_lines; $i++ )
	{
	  print OUT $tail[$i];
	}
      close OUT;
    }
    print "\n";
    unlink "$ARGV[0]/$file.mol2";
}
