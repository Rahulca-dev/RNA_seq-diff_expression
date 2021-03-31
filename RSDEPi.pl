#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
my $gtf_File = "";
my $genome_File = "";
my $dir = "";
my $out_dir = "";


my $result = GetOptions(
	"gtf=s" => \$gtf_File,
	"genome=s" => \$genome_File,
	"o=s" => \$out_dir,
	"raw_data=s" => \$dir,	
);



opendir (my $dh, $dir) || (die "Could not open '$dir' for reading '$!'\n");
my @things = grep {/\.fastq$/} readdir $dh;
my $len = scalar (@things);
print "$len\n";
my %file_hash=();
if($len == 2)
{
	print "Running Single End Procedure\n";
	my $fastq1 = "";
	my $fastq2 = "";
	foreach my $thing (@things) 
	{	
		my $fn = basename($dir."/".$thing);
		print "$fn\t$dir/$thing\n";
		$file_hash{$fn} = $dir."/".$thing;
		if($fastq1 eq "")
		{
			$fastq1 = $thing;
		}
		else
		{
			if($fastq2 eq "")
			{
				$fastq2 = $thing;
			}
		}
	}

	print "RAW DATA FILES:$fastq1\t$fastq2\n";
	my @arr1= split("_",$fastq1);
	my @arr2= split("_",$fastq2);
	if($arr1[0] ne $arr2[0])
	{	
		my$flag1=0;
		my$flag2=0;
	 	if($fastq1=~m/(\w+)\_1/)
		{
			$flag1=1;
		}
		if($fastq2=~m/(\w+)\_1/)
		{
			$flag2=1;
		}
		else
		{	
			print "Not accepted";
		}
		if($flag1==1 && $flag2==1)
		{
			validatefastqcfile($dir."/".$fastq1,$dir."/".$fastq2);

			fastqc($dir."/".$fastq1,$dir."/".$fastq2);
			my @newarr1 = split (/\./,$fastq1);
			my @newarr2 = split (/\./,$fastq2);
			my $file1 = "$out_dir/".$newarr1[0]."_fastqc/fastqc_data.txt";
			my $file2 = "$out_dir/".$newarr2[0]."_fastqc/fastqc_data.txt";
			getOverRepresentedSeq($file1,$file2);
			my $direc1 = dirname($file1);
			my $direc2 = dirname($file2);
			if(-e "$direc1/overrepresented.fasta")
			{
				cutadapt("$direc1/overrepresented.fasta",$direc1,$dir."/".$fastq1,$newarr1[0]);
			}
			if(-e "$direc2/overrepresented.fasta")
			{
				cutadapt("$direc2/overrepresented.fasta",$direc2,$dir."/".$fastq2,$newarr2[0]);
			}
			if(-e "$direc1/trimmed_$newarr1[0].fastq")
			{
				my $trimfastq1 ="$direc1/trimmed_$newarr1[0].fastq";
				fastqc($trimfastq1);
			}
			if(-e "$direc2/trimmed_$newarr2[0].fastq")
			{
				my $trimfastq2 ="$direc2/trimmed_$newarr2[0].fastq";
				fastqc($trimfastq2);
			}
			my $bowtie_dir = $out_dir."/bowtie2DB";
			mkdir $bowtie_dir;	
			Bowtie($genome_File, $bowtie_dir."/genome");
			my $tophat_dir = $out_dir."/tophat_outdir";
			mkdir $tophat_dir;
			my $tophat1 = $tophat_dir."/"."tophat_".$newarr1[0];
			mkdir $tophat1;
			tophat($tophat1,$gtf_File,$bowtie_dir."/genome",$dir."/".$fastq1);
			my $tophat2 = $tophat_dir."/"."tophat_".$newarr2[0];
			mkdir $tophat2;
			tophat($tophat2,$gtf_File,$bowtie_dir."/genome",$dir."/".$fastq2);
			my $cufflinks_dir = $out_dir."/cufflinks_dir";
			mkdir $cufflinks_dir;
			my $cufflinks1 = $cufflinks_dir."/"."cufflinks_".$newarr1[0];
			mkdir $cufflinks1;
			cufflinks($cufflinks1,$gtf_File,$tophat1."/accepted_hits.bam");
			#$file_hash{$fastq1} = $tophat1."/accepted_hits.bam";
			my $cufflinks2 = $cufflinks_dir."/"."cufflinks_".$newarr2[0];
			mkdir $cufflinks2;
			cufflinks($cufflinks2,$gtf_File,$tophat2."/accepted_hits.bam");
			#$file_hash{$fastq2} = $tophat2."/accepted_hits.bam";
			open(my $fw1,">",$out_dir."/assemblies.txt");
			print $fw1 "$cufflinks1/transcripts.gtf\n$cufflinks2/transcripts.gtf";
			my $cuffmerge_dir = $out_dir."/cuffmerge_dir";
			mkdir $cuffmerge_dir;
			cuffmerge($cuffmerge_dir,$gtf_File,$genome_File,$out_dir."/assemblies.txt");
			my $cuffdiff_dir = $out_dir."/cuffdiff_dir";
			mkdir $cuffdiff_dir;
			$file_hash{$fastq1} = $tophat1."/accepted_hits.bam";
			$file_hash{$fastq2} = $tophat2."/accepted_hits.bam";
			open(my $condition_file , "$dir/condition.txt");
			my $c_group = "";
			my $t_group = "";
			my $c_flag = 0;
			my $t_flag = 0;

			while(my $line = <$condition_file>)
			{
				chomp($line);
				my @condition_line1 = split ("\t",$line);
				
				if($condition_line1[1] eq "c")
				{
					if($t_flag == 0)
					{
					 	$c_flag = 1;
					}
					$c_group = $c_group.$file_hash{$condition_line1[0]}.",";		
				}
				if($condition_line1[1] eq "t")
				{
					if($c_flag == 0)
					{
						$t_flag = 1;
					}
					$t_group = $t_group.$file_hash{$condition_line1[0]}.",";
				}
			}

			if($c_flag == 1)
			{
				chop($c_group);chop($t_group);
				print "Control Group : $c_group $t_group\n";
				cuffdiff($cuffdiff_dir,"C,T",$gtf_File,$c_group,$t_group);
			}
			if($t_flag == 1)
			{
				chop($c_group);chop($t_group);	
				print "Treated Group : $t_group $c_group\n";
				cuffdiff($cuffdiff_dir,"T,C",$gtf_File,$t_group,$c_group);
			}
		}
	}
	else
	{
		print "Files Rejected";
	}
}
if($len==4)
{
	print "Running Paired End Procedure\n";
	my @test=sort(@things);
	my $fastq1 = "";
	my $fastq2 = "";
	my $fastq3 = "";
	my $fastq4 = "";
	foreach my $thing (@test) 
	{	
		if($fastq1 eq "")
		{
			$fastq1 = $thing;
		}
		else
		{
			if($fastq2 eq "")
			{
				$fastq2 = $thing;
			}
			else
			{
				if($fastq3 eq "")
				{
					$fastq3=$thing;
				}
				else
				{
					if($fastq4 eq "")
					{
						$fastq4=$thing;
					}
				}				
			}
		}
	}
	print "$fastq1\n$fastq2\n$fastq3\n$fastq4\n";
	my @arr1= split("_",$fastq1);
	my @arr2= split("_",$fastq2);
	my @arr3= split("_",$fastq3);
	my @arr4= split("_",$fastq4);
	print "$arr1[0]\t$arr2[0]\t$arr3[0]\t$arr4[0]\n";
	my $flagg1=0;
	my $flagg2=0;
	my $flagg3=0;
	my $flagg4=0;
	if($arr1[0] eq $arr2[0]||$arr3[0] eq $arr4[0])
	{	
		
	 	if($fastq1=~m/(\w+)\_1/) 
		{
			$flagg1=1;
		}
		if($fastq2=~m/(\w+)\_2/)
		{
			$flagg2=1;
		}	
		if($fastq3=~m/(\w+)\_1/)
		{
			$flagg3=1;	
		}
		if($fastq2=~m/(\w+)\_2/)
		{
			$flagg4=1;
		}
		if($flagg1==1 && $flagg2==1 && $flagg3==1 && $flagg4==1)
		{
			validatefastqcfile($dir."/".$fastq1,$dir."/".$fastq2,$dir."/".$fastq3,$dir."/".$fastq4);
			fastqc($dir."/".$fastq1,$dir."/".$fastq2,$fastq2,$dir."/".$fastq3,$dir."/".$fastq4);
			my @newarr1 = split (/\./,$fastq1);
			my @newarr2 = split (/\./,$fastq2);
			my @newarr3 = split (/\./,$fastq3);
			my @newarr4 = split (/\./,$fastq4);
			my $file1 = "$out_dir/".$newarr1[0]."_fastqc/fastqc_data.txt";
			my $file2 = "$out_dir/".$newarr2[0]."_fastqc/fastqc_data.txt";
			my $file3 = "$out_dir/".$newarr3[0]."_fastqc/fastqc_data.txt";
			my $file4 = "$out_dir/".$newarr4[0]."_fastqc/fastqc_data.txt";
			getOverRepresentedSeq($file1,$file2,$file3,$file4);
			my $direc1 = dirname($file1);
			my $direc2 = dirname($file2);
			my $direc3 = dirname($file3);
			my $direc4 = dirname($file4);
			if(-e "$direc1/overrepresented.fasta")
			{
				cutadapt("$direc1/overrepresented.fasta",$direc1,$dir."/".$fastq1,$newarr1[0]);
			}
			if(-e "$direc2/overrepresented.fasta")
			{
				cutadapt("$direc2/overrepresented.fasta",$direc2,$dir."/".$fastq2,$newarr2[0]);
			}
			if(-e "$direc3/overrepresented.fasta")
			{
				cutadapt("$direc3/overrepresented.fasta",$direc3,$dir."/".$fastq3,$newarr3[0]);
			}
			if(-e "$direc4/overrepresented.fasta")
			{
				cutadapt("$direc4/overrepresented.fasta",$direc4,$dir."/".$fastq4,$newarr4[0]);
			}
			if(-e "$direc1/trimmed_$newarr1[0].fastq")
			{
				my $trimfastq1 ="$direc1/trimmed_$newarr1[0].fastq";
				fastqc($trimfastq1);
			}
			if(-e "$direc2/trimmed_$newarr2[0].fastq")
			{
				my $trimfastq2 ="$direc2/trimmed_$newarr2[0].fastq";
				fastqc($trimfastq2);
			}
			if(-e "$direc3/trimmed_$newarr3[0].fastq")
			{
				my $trimfastq3 ="$direc3/trimmed_$newarr3[0].fastq";
				fastqc($trimfastq3);
			}
			if(-e "$direc4/trimmed_$newarr4[0].fastq")
			{
				my $trimfastq4 ="$direc4/trimmed_$newarr4[0].fastq";
				fastqc($trimfastq4);
			}
			my $bowtie_dir = $out_dir."/bowtie2DB";
			mkdir $bowtie_dir;	
			Bowtie($genome_File, $bowtie_dir."/genome");
			my $tophat_dir = $out_dir."/tophat_outdir";
			mkdir $tophat_dir;
			my $tophat1 = $tophat_dir."/"."tophat_".$newarr1[0];
			mkdir $tophat1;
			tophat($tophat1,$gtf_File,$bowtie_dir."/genome",$dir."/".$fastq1);
			my $tophat2 = $tophat_dir."/"."tophat_".$newarr2[0];
			mkdir $tophat2;
			tophat($tophat2,$gtf_File,$bowtie_dir."/genome",$dir."/".$fastq2);
			my $tophat3 = $tophat_dir."/"."tophat_".$newarr3[0];
			mkdir $tophat3;
			tophat($tophat3,$gtf_File,$bowtie_dir."/genome",$dir."/".$fastq3);
			my $tophat4 = $tophat_dir."/"."tophat_".$newarr4[0];
			mkdir $tophat4;
			tophat($tophat4,$gtf_File,$bowtie_dir."/genome",$dir."/".$fastq4);
			my $cufflinks_dir = $out_dir."/cufflinks_dir";
			mkdir $cufflinks_dir;
			my $cufflinks1 = $cufflinks_dir."/"."cufflinks_".$newarr1[0];
			mkdir $cufflinks1;
			cufflinks($cufflinks1,$gtf_File,$tophat1."/accepted_hits.bam");
			my $cufflinks2 = $cufflinks_dir."/"."cufflinks_".$newarr2[0];
			mkdir $cufflinks2;
			cufflinks($cufflinks2,$gtf_File,$tophat2."/accepted_hits.bam");
			my $cufflinks3 = $cufflinks_dir."/"."cufflinks_".$newarr3[0];
			mkdir $cufflinks3;
			cufflinks($cufflinks3,$gtf_File,$tophat3."/accepted_hits.bam");
			my $cufflinks4 = $cufflinks_dir."/"."cufflinks_".$newarr4[0];
			mkdir $cufflinks4;
			cufflinks($cufflinks4,$gtf_File,$tophat4."/accepted_hits.bam");
			open(my $fw1,">",$out_dir."/assemblies.txt");
			print $fw1 "$cufflinks1/transcripts.gtf\n$cufflinks2/transcripts.gtf\n$cufflinks3/transcripts.gtf\n$cufflinks4/transcripts.gtf";
			my $cuffmerge_dir = $out_dir."/cuffmerge_dir";
			mkdir $cuffmerge_dir;
			cuffmerge($cuffmerge_dir,$gtf_File,$genome_File,$out_dir."/assemblies.txt");
			my $cuffdiff_dir = $out_dir."/cuffdiff_dir";
			mkdir $cuffdiff_dir;	
			$file_hash{$fastq1} = $tophat1."/accepted_hits.bam";
			$file_hash{$fastq2} = $tophat2."/accepted_hits.bam";
			$file_hash{$fastq3} = $tophat3."/accepted_hits.bam";
			$file_hash{$fastq4} = $tophat4."/accepted_hits.bam";
			open(my $condition_file , "$dir/condition.txt");
			my $c_group = "";
			my $t_group = "";
			my $c_flag = 0;
			my $t_flag = 0;

			while(my $line = <$condition_file>)
			{
				chomp($line);
				my @condition_line1 = split ("\t",$line);
				
				if($condition_line1[1] eq "c")
				{
					if($t_flag == 0)
					{
					 	$c_flag = 1;
					}
					$c_group = $c_group.$file_hash{$condition_line1[0]}.",";		
				}
				if($condition_line1[1] eq "t")
				{
					if($c_flag == 0)
					{
						$t_flag = 1;
					}
					$t_group = $t_group.$file_hash{$condition_line1[0]}.",";
				}
			}

			if($c_flag == 1)
			{
				chop($c_group);chop($t_group);
				print "Control Group : $c_group $t_group\n";
				cuffdiff($cuffdiff_dir,"C,T",$gtf_File,$c_group,$t_group);
			}
			if($t_flag == 1)
			{
				chop($c_group);chop($t_group);	
				print "Treated Group : $t_group $c_group\n";
				cuffdiff($cuffdiff_dir,"T,C",$gtf_File,$t_group,$c_group);
			}
		}
	}
	else
	{
		print "Files Rejected";
	}
}

sub getOverRepresentedSeq
{
		
	foreach my $file (@_)
	{
		my $direc = dirname($file);
		my $flagg=0;
		my $header=1;
		open (my $fh,$file)|| die ("\nUnable to open file...\n");
		open (my $fw,'>',"$direc/overrepresented.fasta");
		while (my $linefile=<$fh>)
		{
				if($linefile=~m/\>\>END\_MODULE/)
				{
					$flagg=0;
				}
				if($flagg==1)
				{
					if($linefile =~ m/^#/)	
					{
						
					}
					else
					{
						my @fastarr=split("\t",$linefile);
						print $fw ">Header_$header\n$fastarr[0]\n";
						$header++;
					}
				}				

				if($linefile=~m/>>Overrepresented\ssequences\tfail/)
				{	
						
					$flagg=1;
				}				
		}
	}
}

sub validatefastqcfile
{
	foreach my $fh (@_)
	{
		my $line_position = 0;
		my $count = 0; 
		open (my $open,$fh)|| die ("File not found");
		while (my $line= <$open>)
		{
			 $line_position++;
			if($line=~m/@/)
			{
				$line_position = 1;	
			}
			if($line_position==2)
			{
				if($line=~m/(\w+)/)
				{
				}	
			}
			if($line_position==3)
			{	
				if($line=~m/\+/)
				{
				}
			}
			if($line_position==4)
			{
				if($line=~m/(\w+)/)
				{	
				}	
			}			
		}
				
	}
	print "All Files Validated\n";
}
sub fastqc
{
	print "Running Fastqc.....\n";
	system("fastqc --extract -q @_ -o $out_dir");
	print "Fastqc Successful\n";
}
sub cutadapt
{
	print "Running Cutadapt.....\n";	
	system("cutadapt -b file:$_[0] -o $_[1]/trimmed_$_[3].fastq $_[2] --quiet");
	print "Cutadapt Successful\n";
}
sub Bowtie
{
	print "Running Bowtie2.....\n";
	system("bowtie2-build $_[0] $_[1] --quiet");
	print "Bowtie2 Sucessful\n";
}
sub tophat
{
	system("tophat2 -o $_[0] -G $_[1] $_[2] $_[3]");
}

sub cufflinks
{
	print "Running Cufflinks...\n";
	system("cufflinks -o $_[0] -G $_[1] $_[2] -q");
	print "cufflinks Sucessful\n";
}
sub cuffmerge
{
	print "Running Cuffmerge...\n";
	system("cuffmerge -o $_[0] -g $_[1] -s $_[2] $_[3]");
	print "cuffmerge Sucessful\n";
}
sub cuffdiff
{
	print "Running Cuffdiff...\n";
	system("cuffdiff -o $_[0] -L $_[1] -u $_[2] $_[3] $_[4] -q");
	print "cuffdiff Sucessful\n";
}


