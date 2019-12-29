#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use lib $Bin;

my ($input , $help);
GetOptions(
	"i|input=s"	=>	\$input,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my ($bsample , $bchr , $bloc , $bref , $balt , $bdep , $brate , $btype , $binfo) = ('') x 9;
my $id = '';
my $n = 1;
my $pos = '';
open IN , "$input";
<IN>;
#my @sorted_infos;
my @aa;
while (<IN>){
	chomp;
	my @tmp = split("\t",$_);
	push @aa,[@tmp];
}
my @sorted_infos = sort {$a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]} @aa;

for (@sorted_infos){
	chomp;
	#my @F = split /\t/ , $_;
	my @F = @$_;
	#somatic_C1803125_C1803125-T     chr1    27056326        27056326        C       T 
	my ($sample , $chr , $loc , $ref , $alt , $dep , $dalt) = @F[0,1,2,4,5,13,14];
	my $type;
	if (length($ref) == 1 and length($alt) == 1 and  !($ref=~/-/) and !($alt=~/-/)){
		$type = 'SNP';
	}else{
		$type = 'INDEL';
	}
	my @ti = split /:/ , $F[9];
	#$dep += $dalt;
	
	my $brefl = length($bref);
	my $bloce = $bloc+$brefl-1 if $bloc;
	my $rate = $dalt/$dep;
	my @rate = sort {$a<=>$b} ($rate , $brate) if $brate;
	if ($bsample ne $sample or $bchr ne $chr
			or ($loc-$bloce > 5 or $loc-$bloce < 0) ){
		if ($pos){
			print "$pos\n";
			$pos = '';
		}
	}else{
		$pos = join("\t" , ($bsample , $bchr , $bloc , $bref , $balt , $bdep , $brate , $btype , $binfo)) . "\n" unless $pos;
		$pos .= join("\t" , ($sample , $chr , $loc , $ref , $alt , $dep , $rate , $type , $F[10])) . "\n";
	}
	($bsample , $bchr , $bloc , $bref , $balt , $bdep , $brate , $btype , $binfo) = ($sample , $chr , $loc , $ref , $alt , $dep , $rate , $type , $F[10]);
}
close IN;
if ($pos){
	print "$pos\n";
}


sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: combinePosInVcf.pl
#
#        USAGE: ./combinePosInVcf.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wang Ruiru (wangrr), wangruiru\@berrygenomics.com
# ORGANIZATION: Berry Genomics
#      VERSION: 1.0
#      CREATED: 07/23/18 15:35:22
#     REVISION: ---
#===============================================================================
EOF!
}



