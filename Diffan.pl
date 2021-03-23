#!/usr/bin/perl

use warnings;
use Statistics::R;
use File::Basename;
use Getopt::Long;
my $diff_file ="";
my $out_dir = "";


my $result = GetOptions(
	"diff=s" => \$diff_File,
	"o=s" => \$out_dir	
);

my $R = Statistics::R->new() ;
  
$R-> startR;
$R -> file = $diff_file
$R -> send('df <- read.csv(file, sep="	", fill = TRUE, header = TRUE);');
$R -> send('df1 <- subset(df, p_value < 0.005)');
$R -> send('df2 <- subset(df, p_value < 0.005 && log2.fold_change. > 0.8 )');
$R -> send('df3 <- subset(df, p_value < 0.005 && log2.fold_change. < -0.8 )');
$R -> send('write.csv(df1, "/mnt/e/Bionfo/New folder/diff_gene.txt");');
$R -> send('write.csv(df2, "/mnt/e/Bionfo/New folder/upreg.txt");');
$R -> send('write.csv(df3, "/mnt/e/Bionfo/New folder/downreg.txt");');


