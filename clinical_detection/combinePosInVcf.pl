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

my ($input , $help , $bam);
GetOptions(
	"i|input=s"	=>	\$input,
	"b|bam=s"	=>	\$bam,
	"help"	=>	\$help,
);
my $samtools = "/share/public/software/samtools-1.4/samtools";
$samtools = 'samtools' unless -s $samtools;

if ($help or ! $input){
	&help;
	exit;
}

$/ = '>';
my %fa;
open IN , "/share/work3/capsmart/pipeline/capSMART/CAPcSMART/capSMART/cSMART170503/reference/hg19.fasta";
<IN>;
while (<IN>){
	chomp;
	my $seq = $_;
	s/^(\S+).*\n//;
	my $id = $1;
	s/\n//g;
	$fa{$id} = $_;
}
close IN;
$/ = "\n";

my ($bchr , $bref , $balt , $btype , $binfo) = ('') x 5;
my ($bloc , $baltDep , $brate) = (0) x 3;
my $id = '';
my $n = 1;
open IN , "$input";
while (<IN>){
	if (/^#/){
		print $_;
		next;
	}
	chomp;
	my @F = split /\t/ , $_;
	my ($chr , $loc , $ref , $alt) = @F[0,1,3,4];
	my $type;
	if (length($ref) == 1 and length($alt) == 1){
		$type = 'SNP';
	}else{
		$type = 'INDEL';
		if (length($ref) > length($alt) and length($alt) == 1){
			$type = 'del';
		}elsif (length($ref) < length($alt) and length($ref) == 1 and substr($ref , 0 , 1) eq substr($alt , 0 , 1)){
			$type = 'ins';
		}
	}
	my @ti = split /:/ , $F[9];
	my ($dep , $dalt) = split /,/ , $ti[1];
	my $minalt = (sort {$a<=>$b} ($dalt , $baltDep))[0];
	#$dep += $dalt;
	
	my $brefl = length($bref);
	my $bloce = $bloc+$brefl-1 if $bloc;
	my $rate = $dalt/$dep;
	my @rate = sort {$a<=>$b} ($rate , $brate) if $brate;
	if ($bchr ne $F[0]
		or abs($rate-$brate) > 0.05 
		or $rate[1]/$rate[0] > 1.3
		or ($btype eq 'SNP' and $type eq 'SNP')
		or $loc-$bloce > 5
		or $loc-$bloce < 0
		or ($loc-$bloce == 0 and $type ne 'ins')
		#or abs($dalt-$baltDep)/$minalt > 0.1 
		or readInfo([$bchr,$bloc,$bref,$balt] , [$chr,$loc,$ref,$alt])){

		if ($bchr){
			print "$bchr\t$bloc\t$id\t$bref\t$balt\t.\t$binfo\n";
		}
		($bchr , $bloc , $id , $bref , $balt , $baltDep , $brate , $btype , $binfo) = (@F[0..4] , $dalt , $rate , $type , join("\t" , @F[6..$#F]));
		if ($id eq '.'){
			$id = "Mut$n";
			$n++;
		}
	}else{
		my $ib = '';
		if ($loc - $bloce > 1){
			$ib = substr($fa{$bchr} , $bloce , $loc-$bloce-1);
		}elsif ($loc - $bloce == 0 and $type eq 'ins'){
			$ref =~ s/^\w//;
			$alt =~ s/^\w//;
		}
		$bref .= $ib . $ref;
		$balt .= $ib . $alt;
		$btype = 'INDEL';
		$id = "Mutm$n";
		$n++;
	}
}
close IN;
print "$bchr\t$bloc\t$id\t$bref\t$balt\t.\t$binfo\n" if $bchr;

sub dealPos{
	my ($chr , $str , $ref , $alt) = @{$_[0]};

	if ($chr !~ /chr/){
		$chr = "chr$chr";
	}

	while (length($ref) > 0 and length($alt) > 0 and substr($ref , 0 , 1) eq substr($alt , 0 , 1)){
		$ref =~ s/^\w//; 
		$alt =~ s/^\w//;
		$str++ if length($ref) > 0;
	}

	my $rl = length($ref);
	my $al = length($alt);
	my $end = $str + $rl - 1;
	$end = $str if $end < $str;

	my $type = '';
	if ($str == $end and $rl == $al and $rl == 1){
		$type = 'snp';
	}elsif ($rl == 0){
		$type = 'ins';
	}elsif ($al == 0){
		$type = 'del';
	}else{
		$type = 'xxx';
	}

	return ($chr , $str , $end , $ref , $alt , $type , $rl , $al);
}

sub readInfo{
	my ($pos1 , $pos2) = @_;
	my ($chr1 , $str1 , $end1 , $ref1 , $alt1 , $type1 , $rl1 , $al1) = dealPos($pos1);
	my ($chr2 , $str2 , $end2 , $ref2 , $alt2 , $type2 , $rl2 , $al2) = dealPos($pos2);
	if ($type1 eq 'xxx' or $type2 eq 'xxx'){
		return 1
	}

	my $loc = $str1;
	my $type = '';
	my %dd;
	my $depth;
	my ($da1 , $da2 , $da) = (0 , 0 , 0);
	for my $file (($bam)){
		open BAM , "-|" , "$samtools view -F 1024 $file $chr1\:$str1-$end1";
		while (my $sam = <BAM>){
			my @sp = split /\t/ , $sam;
			my $maplen = &maplen($sp[5]);
			unless ($sp[3] <= $str1 and $end1 <= $sp[3]+$maplen-1 and $sp[3] <= $str2 and $end2 <= $sp[3]+$maplen-1){
				next;
			}
			my ($readn , $refn) = (-1 , $sp[3]-1);
			my ($k1 , $k2) = (0 , 0);
			$depth++;
			while ($sp[5] =~ /(\d+)([IDMS])/g){
				my ($d , $w) = ($1 , $2);
				if ($w =~ /S/){
					$readn += $d;
				}elsif ($w =~ /I/){
					my $bseq = substr($sp[9] , $readn+1 , $d);
					#print STDERR "$refn\t$str1\t$bseq,$alt1\n" if $type1 eq 'ins';
					#print STDERR "$refn\t$str2\t$bseq,$alt2\n" if $type2 eq 'ins';
					if ($refn == $str1 and $type1 eq 'ins' and $bseq eq $alt1){
						$da1++;
						$k1 = 1;
					}elsif ($refn == $str2 and $type2 eq 'ins' and $bseq eq $alt2){
						$da2++;
						$k2 = 1;
					}
					$readn += $d;
				}elsif ($w =~ /D/){
					if ($refn+1==$str1 and $d == $rl1 and $type1 eq 'del'){
						$da1++;
						$k1 = 1;
					}elsif ($refn+1==$str2 and $d == $rl2 and $type2 eq 'del'){
						$da2++;
						$k2 = 1;
					}
					$refn += $d;
				}elsif ($w =~ /M/){
					if ($type1 eq 'snp' and $refn+1<=$str1 and $str1<=$refn+$d and substr($sp[9] , $str1-$refn+$readn , $al1) eq $alt1){
						$da1++;
						$k1 = 1;
					}
					if ($type2 eq 'snp' and $refn+1<=$str2 and $str2<=$refn+$d and substr($sp[9] , $str2-$refn+$readn , $al2) eq $alt2){
						$da2++;
						$k2 = 1;
					}
					$refn += $d;
					$readn += $d;
				}else{
					next;
				}
			}
			if ($k1 and $k2){
				$da++;
			}
		}
		close BAM;
	}
	my $dam = (sort {$a<=>$b} ($da1 , $da2))[0];
	print STDERR "$chr1\t$str1\t$end1\t$ref1\t$alt1\t$da1\t$da\n";
	print STDERR "$chr2\t$str2\t$end2\t$ref2\t$alt2\t$da2\t$da\n";
	return 1 if $dam == 0;
	if ($da/$dam < 0.5){
		return 1;
	}else{
		return 0;
	}
}

sub maplen{
	my ($mode) = @_;
	my $len = 0;
	while ($mode =~ /(\d+)([MSHID])/g){
		if ($2 eq 'M' or $2 eq 'D'){
			$len += $1;
		}
	}
	return $len;
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



