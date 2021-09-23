#!/usr/bin/perl -w

#################################################################################
# This is the MERCI (Motif - EmeRging and with Classes - Identification) program.
# Copyright (C) 2010  Celine Vens
# 
# Usage will be shown when script is run without any arguments: "perl MERCI.pl"
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#     
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#     
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#################################################################################    


use strict;
use Time::HiRes;
#use Module::Load;

# GLOBAL VARIABLES
##################

# input arguments
my @arguments = @ARGV;
my $posfile = -1; # fasta file with positive sequences (REQUIRED)
my $minfreq = 1; # minimal frequency threshold (absolute number) for positive sequences (global var)
my $negfile = -1; # fasta file with negative sequences (REQUIRED)
my $maxfreq = 0; # maximal frequency threshold (absolute number) for negative negatives (global var)
my $patternsfile = "motifs"; # output file to store motifs
my $hierarchy = "NONE"; # class hierarchy to use (NONE/KOOLMAN-ROHM/BETTS-RUSSELL/RASMOL/hierarchy_file)
my $maxlength = 10000; # maximal number of elements in the patterns
my $maxlevels = 10000; # maximum number of levels to explore
my $parallel = 0; # 0 when not in parallel mode, root refinement when in parallel mode
my $stringlength = 10000;
my $maxgaplength = 1;
my $maxgaps = 0;
my $k = 10; # top K

# some timers and counters for printing statistics
my $refinetime = 0;
my $amtime = 0;
my $mtime = 0;
my $prunetime = 0;

# some hashes representing the hierarchy
my %string_representation; # the representation of the concepts in the data file
my %minimal_refinements; # list all minimal refinements of this concept that will be generated, ensure TREE
my %DAG_parents; # list all other maximally specific parents (not in the chosen tree) of this concept

# hash representing the top K patterns
my %topk;

# storing the results
my @patterns; # collection of motifs that are found

# print message
my $message;


# PROCESSING ARGUMENTS
######################
print "\nThis is the MERCI (Motif - EmeRging and with Classes - Identification) program.\n\n";

if(@arguments <= 0 ){
	print "To run the program, please define at least the required options:\n\n";
	print "  -p posfile          # the file with positive sequences (REQUIRED)\n";
	print "  -n negfile          # the file with negative sequences (REQUIRED)\n";
	print "  -k topK			 # the number of motifs requested (values: number/ALL, default: 10)\n";
	print "  -fp posfreq         # the minimal frequency (absolute number) for the positive sequences (default: 1)\n";
	print "  -fn negfreq         # the maximal frequency (absolute number) for the negative sequences (default: 0)\n";
	print "  -o outputfile       # file where found motifs will be stored (default: motifs)\n";
	print "  -l maxlength        # the maximal motif length (default: 10000)\n";
	print "  -c classification   # the classification hierarchy to use (values: NONE/KOOLMAN-ROHM/BETTS-RUSSELL/RASMOL/my_classification_file, default: NONE)\n";
	print "  -s sequencelength   # only the first sequencelength positions of the sequences will be considered (default: 10000)\n";
	print "  -g maxnbgaps	     # maximal number of gaps (default: 0)\n";
	print "  -gl maxgaplength    # maximal gap length (default: 1)\n";
	print "  -para parallel      # whether to run only a limited part of the search (values: 0/first_level_refinement, default: 0)\n\n";
	print "See manual for more information on the options.\n\n";
	exit;
}
else
{
	my $option;
	my $value;
	while ($option = shift(@arguments))
	{
		$value = shift(@arguments);
		&process_option($option,$value);
	}
}
&check_options();


sub process_option()
{
	my $option = $_[0];
	my $value = $_[1];
	
	if ($option eq "-p")
	{
		if (not defined $value)
		{
			die "File $value (option -p) not found.\n";
		}
		elsif (not -e $value)
		{
			die "File $value (option -p) not found.\n";
		}
		$posfile = $value;
	}
	elsif ($option eq "-n")
	{
		if (not defined $value)
		{
			die "File $value (option -n) not found.\n";
		}
		elsif (not -e $value)
		{
			die "File $value (option -n) not found.\n";
		}
		$negfile = $value;
	}
	elsif ($option eq "-k")
	{
		$k = $value;
	}
	elsif ($option eq "-fp")
	{
		$minfreq = $value;
	}
	elsif ($option eq "-fn")
	{
		$maxfreq = $value;
	}
	elsif ($option eq "-o")
	{
		$patternsfile = $value;
	}
	elsif ($option eq "-l")
	{
		$maxlength = $value;
	}
	elsif ($option eq "-g")
	{
		$maxgaps = $value;
	}
	elsif ($option eq "-gl")
	{
		$maxgaplength = $value;
	}
	elsif ($option eq "-c")
	{
		if (((((not $value eq "NONE") and (not $value eq "KOOLMAN-ROHM")) and (not $value eq "BETTS-RUSSELL")) and (not $value eq "RASMOL")))
		{
			if (not -e $value)
			{
				die "Classification file $value (option -c) not found. Use an existing file name, or one of NONE/KOOLMAN-ROHM/BETTS-RUSSELL/RASMOL.\n";
			}
		}
		$hierarchy = $value;
	}
	elsif ($option eq "-para")
	{
		$parallel = $value;
	}
	elsif ($option eq "-s")
	{
		$stringlength = $value;
	}
	else
	{
		die "Unknown option $option. Run the program without any options to see its usage.\n";
	}
}

sub check_options()
{
	if ($posfile eq "-1")
	{
		die "No fasta file with positive sequences given. Please set option -p.\n";
	}
	if ($negfile eq "-1")
	{
		die "No fasta file with negative sequences given. Please set option -n.\n";
	}
}


# initialize output file
########################
open(PAT,">$patternsfile");
print PAT "MERCI - Summary of options:\n";
print PAT " -p    $posfile\n";
print PAT " -n    $negfile\n";
print PAT " -k 	  $k\n";
print PAT " -fp   $minfreq\n";
print PAT " -fn   $maxfreq\n";
print PAT " -o    $patternsfile\n";
print PAT " -l    $maxlength\n";
print PAT " -s    $stringlength\n";
print PAT " -c    $hierarchy\n";
print PAT " -g    $maxgaps\n";
print PAT " -gl   $maxgaplength\n";
print PAT " -para $parallel\n\n\n";
print PAT "Motifs:\n";
close(PAT);



# HIERARCHY
###########
print "   * Loading hierarchy $hierarchy\n";

$string_representation{"root"} = 'root';
$string_representation{"gap"} = '.{0,' . $maxgaplength . '}';

$string_representation{"A"} = 'A';
$string_representation{"C"} = 'C';
$string_representation{"D"} = 'D';
$string_representation{"E"} = 'E';
$string_representation{"F"} = 'F';
$string_representation{"G"} = 'G';
$string_representation{"H"} = 'H';
$string_representation{"I"} = 'I';
$string_representation{"K"} = 'K';
$string_representation{"L"} = 'L';
$string_representation{"M"} = 'M';
$string_representation{"N"} = 'N';
$string_representation{"P"} = 'P';
$string_representation{"Q"} = 'Q';
$string_representation{"R"} = 'R';
$string_representation{"S"} = 'S';
$string_representation{"T"} = 'T';
$string_representation{"U"} = 'U';
$string_representation{"V"} = 'V';
$string_representation{"W"} = 'W';
$string_representation{"Y"} = 'Y';

&load_hierarchy($hierarchy);

sub load_hierarchy()
{
	my $hierarchy = $_[0];
	
	if ($hierarchy eq "NONE") # hierarchy without concepts
	{
		$minimal_refinements{"root"} = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"];
	}
	
	elsif ($hierarchy eq "KOOLMAN-ROHM")
	{
		$string_representation{"aliphatic"} = '[AGILV]';
		$string_representation{"sulfur"} = '[CM]';
 		$string_representation{"aromatic"} = '[FYW]';
		$string_representation{"neutral"} = '[STNQ]';
		$string_representation{"acidic"} = '[DE]';
		$string_representation{"basic"} = '[RHK]';
		
		$minimal_refinements{"root"} = ["aliphatic","sulfur","aromatic","neutral","acidic","basic","P"];
		$minimal_refinements{"aliphatic"} = ["A","G","I","L","V"];
		$minimal_refinements{"sulfur"} = ["C","M"];
		$minimal_refinements{"aromatic"} = ["F","Y","W"];
		$minimal_refinements{"neutral"} = ["S","T","N","Q"];
		$minimal_refinements{"acidic"} = ["D","E"];
		$minimal_refinements{"basic"} = ["R","H","K"];
	}
	
	elsif ($hierarchy eq "BETTS-RUSSELL")
	{
		$string_representation{"polar"} = '[HKRDEYWTCSNQ]';
		$string_representation{"charged"} = '[DERHK]';
		$string_representation{"negative"} = '[DE]';
		$string_representation{"positive"} = '[RHK]';
		$string_representation{"small"} = '[AGCSPNDTV]';
		$string_representation{"tiny"} = '[AGCS]';
		$string_representation{"hydrophobic"} = '[HFWYILVMKTAGC]';
		$string_representation{"aromatic"} = '[HFWY]';
		$string_representation{"aliphatic"} = '[ILV]';
		
		$minimal_refinements{"root"} = ["polar","hydrophobic","small"];
		$minimal_refinements{"polar"} = ["charged","Q"];
		$minimal_refinements{"charged"} = ["positive","negative"];
		$minimal_refinements{"positive"} = ["R"];
		$minimal_refinements{"negative"} = ["E"];
		$minimal_refinements{"hydrophobic"} = ["aromatic","aliphatic","M","K"];
		$minimal_refinements{"aromatic"} = ["F","Y","W","H"];
		$minimal_refinements{"aliphatic"} = ["I","L"];
		$minimal_refinements{"small"} = ["tiny","P","N","D","T","V"];
		$minimal_refinements{"tiny"} = ["A","C","G","S"];
		
		$DAG_parents{"A"} = ["hydrophobic"]; #tiny
		$DAG_parents{"C"} = ["hydrophobic","polar"]; #tiny
		$DAG_parents{"D"} = ["negative"]; #small
		$DAG_parents{"G"} = ["hydrophobic"]; #tiny
		$DAG_parents{"H"} = ["positive"]; #aromatic
		$DAG_parents{"K"} = ["positive"]; #hydrophobic
		$DAG_parents{"N"} = ["polar"]; #small
		$DAG_parents{"S"} = ["polar"]; #tiny
		$DAG_parents{"T"} = ["hydrophobic","polar"]; #small
		$DAG_parents{"V"} = ["aliphatic"]; #small
		$DAG_parents{"W"} = ["polar"]; #aromatic
		$DAG_parents{"Y"} = ["polar"]; #aromatic
	}
	
	elsif ($hierarchy eq "RASMOL")
	{
		$string_representation{"charged"} = '[DERHK]';
		$string_representation{"acidic"} = '[DE]';
		$string_representation{"basic"} = '[RHK]';
		$string_representation{"neutral"} = '[ANCQGILMFPSTWYV]';
		$string_representation{"cyclic"} = '[HFWYP]';
		$string_representation{"acyclic"} = '[ARNDCEQGILKMSTV]';
		$string_representation{"aromatic"} = '[HFWY]';
		$string_representation{"aliphatic"} = '[AGILV]';
		$string_representation{"surface"} = '[RNDEQGHKPSTY]';
		$string_representation{"buried"} = '[ACILMFWV]';
		$string_representation{"hydrophobic"} = '[AGILMFPWYV]';
		$string_representation{"polar"} = '[RNDCEQHKST]';
		$string_representation{"small"} = '[AGS]';
		$string_representation{"medium"} = '[NDCPTV]';
		$string_representation{"large"} = '[REQHILKMFWY]';

		$minimal_refinements{"root"} = ["neutral","acyclic","hydrophobic","large","cyclic","surface","polar","medium"];
		$minimal_refinements{"neutral"} = ["buried"];
		$minimal_refinements{"acyclic"} = ["small"];
		$minimal_refinements{"hydrophobic"} = ["aliphatic"];
		$minimal_refinements{"large"} = ["I","L","M"];
		$minimal_refinements{"cyclic"} = ["aromatic"];
		$minimal_refinements{"surface"} = ["G","Y"];
		$minimal_refinements{"polar"} = ["charged","Q","S"];
		$minimal_refinements{"medium"} = ["N","D","C","P","T","V"];
		$minimal_refinements{"aliphatic"} = ["A"];
		$minimal_refinements{"aromatic"} = ["F","W"];
		$minimal_refinements{"charged"} = ["acidic","basic"];
		$minimal_refinements{"basic"} = ["R","H","K"];
		$minimal_refinements{"acidic"} = ["E"];

		$DAG_parents{"A"} = ["buried","small"]; #aliphatic
		$DAG_parents{"C"} = ["buried","polar","acyclic"]; #medium
		$DAG_parents{"D"} = ["acidic"]; #medium
		$DAG_parents{"E"} = ["large"]; #acidic
		$DAG_parents{"F"} = ["buried","hydrophobic"]; #aromatic
		$DAG_parents{"G"} = ["aliphatic","small"]; #surface
		$DAG_parents{"H"} = ["aromatic"]; #basic
		$DAG_parents{"I"} = ["aliphatic","buried","hydrophobic"]; #large
		$DAG_parents{"K"} = ["acyclic"]; #basic
		$DAG_parents{"L"} = ["aliphatic","buried","hydrophobic"]; #large
		$DAG_parents{"M"} = ["acyclic","buried","hydrophobic"]; #large
		$DAG_parents{"N"} = ["neutral","surface","acyclic","polar"]; #medium
		$DAG_parents{"P"} = ["neutral","surface","hydrophobic","cyclic"]; #medium
		$DAG_parents{"Q"} = ["neutral","surface","large","acyclic"]; #polar
		$DAG_parents{"R"} = ["acyclic"]; #basic
		$DAG_parents{"S"} = ["surface","small"]; #polar
		$DAG_parents{"T"} = ["neutral","surface","acyclic","polar"]; #medium
		$DAG_parents{"V"} = ["buried","aliphatic"]; #medium
		$DAG_parents{"W"} = ["buried","hydrophobic"]; #aromatic
		$DAG_parents{"Y"} = ["neutral","hydrophobic","aromatic"]; #surface
		$DAG_parents{"small"} = ["neutral"];
		$DAG_parents{"aliphatic"} = ["neutral","acyclic"];
		$DAG_parents{"aromatic"} = ["large"];
		$DAG_parents{"charged"} = ["surface"];
		$DAG_parents{"basic"} = ["large"];
		$DAG_parents{"acidic"} = ["acyclic"];
	}
	
	else
	{
		&load_user_hierarchy($hierarchy);
	}
}


# MAIN PROGRAM
##############

# read sequences and store them as an array of strings and headers
print "   * Reading data\n";
my ($posset,$posheaders) = &extract_sequences($posfile); # (global var)
my ($negset,$negheaders) = &extract_sequences($negfile); # (global var)

# mine patterns
print "   * Searching for motifs\n";
my $start = [ Time::HiRes::gettimeofday( ) ];
&mine_root();
if (not $k eq "ALL")
{
	&get_topK();
}
my $elapsed = Time::HiRes::tv_interval( $start );
print "\n      " . @patterns . " motifs printed to $patternsfile\n\n\n";

# print time statistics
print "      Total time taken to find motifs: $elapsed seconds\n";
print "       -> Refinement: $refinetime seconds\n";
print "       -> Pruning: $prunetime seconds\n";
print "       -> Checking occurrence in positives: $amtime seconds\n";
print "       -> Checking occurrence in negatives: $mtime seconds\n\n";

# print time statistics to output file
# open(PAT,">>$patternsfile");
# print PAT "\n\n===================================================\n";
# print PAT "       Total time taken to find motifs: $elapsed seconds\n";
# print PAT "       -> Refinement: $refinetime seconds\n";
# print PAT "       -> Pruning: $prunetime seconds\n";
# print PAT "       -> Checking occurrence in positives: $amtime seconds\n";
# print PAT "       -> Checking occurrence in negatives: $mtime seconds\n";
# close(PAT);

# find the locations of the patterns
print "   * Looking for occurrences of found motifs\n";
&check_locations(@patterns);


# SUBROUTINES
#############

sub mine_root()
{
	my $rootstring = "root";
	my $candidate = [$rootstring];
	
	$minimal_refinements{"rootall"} = $minimal_refinements{"root"};
	if (not $parallel eq "0") # if parallel option: change $minimal_refinements{"root"}
	{
		$minimal_refinements{"root"} = [ $parallel ];
	}
	
	my $refinements = &refine_candidate($candidate);
	my @examples = (0 .. @$posset-1);
	
	foreach my $ref (@$refinements)
	{
		my %amhash = ();
		&mine($ref,\@examples,1,\%amhash,0);
	}
}

sub mine()
{
	my $candidate = $_[0]; # candidate to be checked, and possibly be refined
	my $examples = $_[1]; # the list of positive sequences to check
	my $level = $_[2]; # refinement level
	my $amhash = $_[3]; # hash to collect patterns that are frequent in the positives
	my $parentmonotone = $_[4]; # if the candidate is infrequent in the negatives, as defined by its parents
	
	my $candstring = &pattern_to_string($candidate);
	# print "checking $candstring\n";

	my $start;
	my $monotone = 0;
	
	$start = [ Time::HiRes::gettimeofday( ) ];
	my $examplestocheck = &prune($candidate,$amhash,$examples);
	$prunetime += Time::HiRes::tv_interval( $start );
	
	if (@$examplestocheck >= $minfreq)
	{
		$start = [ Time::HiRes::gettimeofday( ) ];
		my $newexamples = &check_antimonotone($candidate,$examplestocheck,$amhash);
		my $nbnewexamples = @$newexamples;
		$amtime += Time::HiRes::tv_interval( $start );
		
		if ($nbnewexamples >= $minfreq) # frequent in the positives
		{
			${$amhash}{$candstring} = [ @$newexamples ];

			my $negatives_covered = -1;
			if ($parentmonotone < 1)
			{
				$start = [ Time::HiRes::gettimeofday( ) ];
				$negatives_covered = &check_monotone($candidate);
				$mtime += Time::HiRes::tv_interval( $start );
			}
			if ($negatives_covered <= $maxfreq) # infrequent in the negatives
			{
				$monotone = 1;
				&add_pattern($candidate,$nbnewexamples);
				&print_pattern($candidate,$nbnewexamples,$negatives_covered);
			}
			
			$start = [ Time::HiRes::gettimeofday( ) ];
			my ($addrefinements,$specrefinements) = &refine_candidate($candidate,$amhash);
			$refinetime += Time::HiRes::tv_interval( $start );
			
			$level++;
			if ($level <= $maxlevels)
			{
				foreach my $ref (@$specrefinements)
				{
					&mine($ref,$newexamples,$level,$amhash,$monotone);
				}
				my %amhash2 = (); # we can start a new hash to save on memory
				foreach my $ref (@$addrefinements)
				{
					&mine($ref,$newexamples,$level,\%amhash2,$monotone);
				}
			}
		}
	}
}


sub add_pattern()
{
	my $candidate = $_[0]; # the new valid pattern found
	my $posfreq = $_[1]; # the frequency in the positives of the new candidate
	# opti: telkens als een candidate in de hash komt: remove generals
	if ($k eq "ALL") # report all motifs with frequency in positives >= $minfreq and frequency in negatives >= $maxfreq
	{
		push(@patterns,$candidate);
		&print_pattern_to_file($candidate);
	}
	else
	{
		if (exists($topk{$posfreq}))
		{
			push @{ $topk{$posfreq} }, $candidate;
		}
		else
		{
			$topk{$posfreq}[0] = $candidate;
		}
		# check if minfreq has to be augmented
		my @frequencies = sort {$b <=> $a} keys %topk;
		my $count = 0;
		my $i = 0;
		my $nbfrequencies = @frequencies;
		while (($count < $k) and ($i<$nbfrequencies))
		{
			$count += @{ $topk{$frequencies[$i]} };
			$i++;
		}
		if ($count >= $k)
		{
			$i--;
			$minfreq = $frequencies[$i];
		}
	}
}

sub get_topK()
{
	my @frequencies = sort {$b <=> $a} keys %topk;
	my $i = 0;	
	my $nbfrequencies = @frequencies;
	while (($i<$nbfrequencies) and ($frequencies[$i] >= $minfreq))
	{
		foreach my $candidate (@{ $topk{$frequencies[$i]} })
		{
			push(@patterns,$candidate);
			&print_pattern_to_file($candidate);
		}
		$i++;
	}
}


# check whether all minimal DAG parents of the last element are frequent in the positives
sub prune()
{
	my $patt = $_[0]; # ref to array of strings
	my $amhash = $_[1];
	my $parentexamples = $_[2];
	my @pattern = @{$patt}; # array of strings
	if ((@pattern > 1) and (exists $DAG_parents{$pattern[@pattern-1]}))
	{
		my @examplestocheck = ();
		my @examplecounters = ();
		for (my $i=0;$i<@$posset;$i++)
		{
			$examplecounters[$i] = 0;
		}
		
		foreach my $ind (@$parentexamples)
		{
			$examplecounters[$ind]++;
		}
		
		foreach my $gen (@{ $DAG_parents{$pattern[@pattern-1]} })
		{
			my @parent_pattern = @pattern;
			$parent_pattern[@pattern-1] = $gen;
			my $parent_string = &pattern_to_string(\@parent_pattern);
			if (not exists $$amhash{$parent_string})
			{
				return [];
			}
			else
			{
				foreach my $ind (@{$$amhash{$parent_string}})
				{
						$examplecounters[$ind]++;
				}
			}
		}
		
 		# pruning is done, now compute the examples to check
		my $nbparentschecked = @{ $DAG_parents{$pattern[@pattern-1]} } + 1;
		for (my $i=0; $i<@examplecounters; $i++)
		{
			if ($examplecounters[$i] == $nbparentschecked)
			{
					push(@examplestocheck,$i);
			}
		}
		return \@examplestocheck;
	}
	else
	{
		return $parentexamples;
	}
}


# print found motif both on standard output and in motif file
sub print_pattern()
{
	my $p = $_[0];
	my $posfreq = $_[1];
	my $negfreq = $_[2];
	print "     * motif: ";
	foreach my $pp (@$p)
	{
		print " $pp";
	}
	print " ($posfreq,$negfreq) (minfreq = $minfreq)\n";
}

sub print_pattern_to_file()
{
	my $p = $_[0];
	open(PAT,">>$patternsfile");
	foreach my $pp (@$p)
	{
		print PAT " $pp";
	}
	print PAT "\n";
	close(PAT);
}


# read the input file and return the sequences as an array of strings
sub extract_sequences()
{
	my $file = $_[0];
	open(IN,"$file") or die "file $file not found\n";
	my @lines = <IN>;
	close(IN);
	my @sequenceset;
	my @sequenceheaders;
	my $sequence = "";
	foreach my $l (@lines)
	{
		if ($l =~ /^>/) # header line
		{
			if (not $sequence eq "") # if not first line
			{
				my $substring = substr($sequence,0,$stringlength);
				push (@sequenceset, $substring);
			}
			$sequence = "";
			push (@sequenceheaders,$l);
		}
		else # sequence line, no header line
		{
			chomp($l);
			$l =~ s/^\s+//; # remove whitespace in front and after sequence
			$l =~ s/\s+$//;
			$sequence = $sequence . $l;
		}
	}
	my $substring = substr($sequence,0,$stringlength);
	push (@sequenceset, $substring); # also add last line
	print "       input file $file read -> " . @sequenceset . " sequences found\n";
	return (\@sequenceset,\@sequenceheaders);
}

# refine one candidate pattern
sub refine_candidate()
{
	my $cand = $_[0]; # reference to array of strings
	my $amhash = $_[1];
	
	my $addrefinements;
	my $specrefinements;
	
	$specrefinements = &specialize($cand,$amhash);
	if (@{$cand} < $maxlength)
	{
		$addrefinements = &add($cand);
	}

	return ($addrefinements,$specrefinements);
}


# add one element of first level to end of pattern
# if pattern is "root" (beginning of mining process), don't add anything (otherwise we get duplicated patterns)
sub add()
{
	# input: ref to array strings
	my $patt = $_[0]; # ref to array of strings
	my @pattern = @{$patt}; # array of strings
	my @refinements = ();
	if (not $pattern[0] eq "root")
	{
		my @firstlevelelements = @{ $minimal_refinements{"rootall"} };
		foreach my $fle (@firstlevelelements)
		{
			# without gap
			my @newpattern = @{$patt};
			push(@newpattern,$fle);
			push(@refinements,[ @newpattern ]);
			# with gap
			if ($maxgaps > 0)
			{
 				if ((&number_gaps($patt) < $maxgaps) and (@pattern < $maxlength-1))
# 				if (@pattern < $maxlength-1)
				{
					@newpattern = @{$patt};
					push(@newpattern,"gap");
					push(@newpattern,$fle);
					push(@refinements,[ @newpattern ]);
				}
			}
		}
	}
	return \@refinements;
}

sub contains_gap()
{
	my $patt = $_[0]; # ref to array of strings
	foreach my $el (@{$patt})
	{
		if ($el eq "gap")
		{
			return 1;
		}
	}
	return 0;
}

sub number_gaps()
{
	my $patt = $_[0]; # ref to array of strings
	my $nbgaps = 0;
	foreach my $el (@{$patt})
	{
		if ($el eq "gap")
		{
			$nbgaps++;
		}
	}
	return $nbgaps;
}

# minimally specialize last element of pattern
sub specialize()
{
	# input: ref to array of strings
	my $patt = $_[0]; # ref to array of strings
	my @pattern = @{$patt}; # array of strings
	my @refinements = ();
	my $lastel = pop(@pattern);
	my @minimalrefs = ();
	if (defined $minimal_refinements{$lastel})
	{
		@minimalrefs = @{ $minimal_refinements{$lastel} };
		foreach my $ref (@minimalrefs)
		{
			my @newpattern = @pattern;
			push(@newpattern,$ref);
			push(@refinements,[ @newpattern ]);
		}
	}
	return \@refinements;
}


# check the positive sequences for occurrence of the candidate (stop counting if minfreq can not be obtained)
sub check_antimonotone()
{
	my $cand = $_[0]; # the candidate to be checked
	my $examples = $_[1]; # the sequences to check coverage (array of indices)
	my $amhash = $_[2]; # hash to collect am alerts
	
	my $count = 0; # the coverage
	my $candstring = &pattern_to_string($cand); # the candidate in string format, with the concepts replaced by a regular expression
	my @newexampleset = (); # the new sequence ids to be checked by the refinements of this candidate in the next level
	
	my $nbtocheck = @{$examples};
	while (($nbtocheck >= $minfreq - $count) and ($nbtocheck > 0)) # optimization: stop checking antimonotonicity if e.g. still 20 hits left to reach minfreq and only 19 sequences left to check
	{		
		my $sequence = $$posset[${$examples}[$nbtocheck-1]];
		if ($sequence =~ /$candstring/) # if candidate pattern found in sequence, increase count, and add id to the new sequence list
		{
			$count++;
			push(@newexampleset,${$examples}[$nbtocheck-1]);
		}
		$nbtocheck--;
	}
	return \@newexampleset;
}


# check the number of negative sequences containing the candidate (stop counting when maxfreq is obtained)
sub check_monotone()
{
	my $cand = $_[0]; # the candidate to be checked
	
	my $candstring = &pattern_to_string($cand); # the candidate in string format, with the concepts replaced by a regular expression
	my $count = 0; # optimization: from the moment we get above $maxfreq, stop counting (fails)
	my $index = 0; # index of negative sequence being checked
	
	while (($index<@{$negset}) and ($count<=$maxfreq)) # stop searching if maxfreq obtained 
	{		
		my $neg = $$negset[$index];
		if ($neg=~/$candstring/) # if candidate pattern found in sequence, increase count
		{
			$count++;
		}
		$index++;
	}
	return $count;
}

# changes a pattern into a string, with the hierarchy concepts replaced by a regular expression
sub pattern_to_string()
{
	my $pattern = $_[0]; # ref to array of strings
	my $patternstring = "";
	foreach my $str (@{$pattern}) # for each string in the array
	{
		if (defined($string_representation{$str}))
		{
			$patternstring = $patternstring . $string_representation{$str};
		}
		else
		{
			die "no string representation for $str\n";
 		}
	}
	return $patternstring;
}

# checking at which locations in the sequences the patterns occur
sub check_locations()
{
	system("perl C:\\Users\\bipas\\MERCI\\MERCI_motif_locator.pl -i $patternsfile -p $posfile -n $negfile -c $hierarchy -gl $maxgaplength -s $stringlength");
}

# load a classification scheme defined by the user (see manual)
sub load_user_hierarchy()
{
	my $file = $_[0];
	open(IN,$file) or die "Can not open hierarchy file $file\n";
	my $line = <IN>;
	while ($line !~ /Definitions/)
	{
		$line = <IN>;
	}
	$line = <IN>;
	while ($line =~ /=/)
	{
		chomp($line);
		$line =~ /(\S+)\s*=\s*(.+)/ or die "Wrong line format in hierarchy file: $line\n";
		my $class = $1;
		my $values = $2;
		$values =~ s/ //g; #remove spaces
		$values =~ s/,//g; #remove commas
		$string_representation{$class} = "\[$values\]";
		$line = <IN>;
	}
	while ($line !~ /Tree structure/)
	{
		$line = <IN>;
	}
	$line = <IN>;
	while (defined($line) and ($line =~ /=/))
	{
		chomp($line);
		$line =~ /(\S+)\s*=\s*(.+)/ or die "Wrong line format in hierarchy file: $line\n";
		my $class = $1;
		my $values = $2;
		$values =~ s/ //g; #remove spaces
		my @allvalues = split(/,/,$values);
		$minimal_refinements{$class} = \@allvalues;
		$line = <IN>;
	}
	while (defined($line) and ($line !~ /DAG parents/))
	{
		$line = <IN>;
	}
	while ($line = <IN>)
	{
		if ($line =~ /=/)
		{
			chomp($line);
			$line =~ /(\S+)\s*=\s*(.+)/ or die "Wrong line format in hierarchy file: $line\n";
			my $class = $1;
			my $values = $2;
			$values =~ s/ //g; #remove spaces
			my @allvalues = split(/,/,$values);
			$DAG_parents{$class} = \@allvalues;
		}
	}
}
