#!perl
# hbind_setup.pl    Volker Schnecke  Wed Nov 17 14:36:37 EST 1999
#
# Usage:  hbind_setup.pl  target template [database]
#
# This script creates the directory structure for a screening experiment.

if ( $#ARGV != 1 && $#ARGV != 2 )
{
    die "Usage: $0 <target> <template> [<database>]\n";
}
 
if ( $ENV{"HBIND_DATA_DIR"} eq "" )
{
    print "ERROR: environment variable HBIND_DATA_DIR not set\n";
    exit;
}

$data_dir = $ENV{"HBIND_DATA_DIR"};
$target = $ARGV[0];
$template = $ARGV[1];
$database = $ARGV[2] if $#ARGV == 2;

print "\nsetting up HBIND screening directories for:\n";
print "  - target protein code: $target\n";
if ( $#ARGV == 1 )
{
    print "  - template specifier:  $template\n\n";
}
else
{
    print "  - template specifier:  $template\n";
    print "  - screening database:  $database\n\n";
}

chdir $data_dir;
if ( ! -e "$target" )
{
    mkdir "$target", 0755;
    print "creating directory $data_dir/$target\n";
}
chdir $target;
if ( ! -e "$template" )
{
    mkdir "$template", 0755;
    print "creating directory $data_dir/$target/$template\n";
}
chdir $template;
if ( ! -e "in" )
{
    mkdir "in", 0755;
    print "creating directory $data_dir/$target/$template/in\n";
}
if ( ! -e "log" )
{
    mkdir "log", 0755;
    print "creating directory $data_dir/$target/$template/log\n";
}
if ( $#ARGV == 2 && ! -e "$database"."_ligands" )
{
    mkdir "$database"."_ligands", 0755;
#    mkdir "$database"."_waters", 0755;
    mkdir "$database"."_targets", 0755;
    print "creating directory $data_dir/$target/$template/$database"."_ligands\n";
#    print "creating directory $data_dir/$target/$template/$database"."_waters\n";
    print "creating directory $data_dir/$target/$template/$database"."_targets\n";
}
print "\n";
