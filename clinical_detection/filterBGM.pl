use strict;
use warnings;
use Statistics::Descriptive;
use Statistics::R;
use Getopt::Long;

my ($in,$out,$info,$help,$bgm_p,$bgm_t,$pvalue);
GetOptions(
  "i|in=s"=>\$in,
  "n|info=s"=>\$info,
  "bp|bgmp=s"=>\$bgm_p,
  "bt|bgmt=s"=>\$bgm_t,
  "o|out=s"=>\$out,
  "p|pvalue=f"=>\$pvalue,
  "h|help"=>\$help,
);

usage() unless $in && $bgm_p && $bgm_t && $info&& $out ;
usage() if $help;

#$pvalue ||= 0.05;
$pvalue ||= 0.05;
sub usage{
	print "usage:$0 <-i in> <-n info> <-o out> <-bp bgmp>  <-bt bgmt>\n";
	exit(1);
}
my %record;
open IN,"$in" or die $!;
my @lines = <IN>;
for (@lines){
	next if /^chrom/;
	next if /strand_filtered$/;
	my @tmp = split("\t",$_);
	if (scalar @tmp >17){
		my @a = @tmp[0..16];
		@tmp = @a;
	}
	my ($chrom,$pos,$ref,$alt) = ($tmp[1],$tmp[2],$tmp[4],$tmp[5]);
	my $key = join("\t",($chrom,$pos,$ref,$alt));
	#print "$key\n";
	$record{$key} = 1;
}
close IN;

my %sampletype;
open INFO,"$info" or die $!;
while (<INFO>){
	my @tmp = split("\t",$_);
	my($subsample, $subtype) = ($tmp[0],$tmp[1]);	
	#print $subsample . " " . $subtype;
	$sampletype{$subsample} = $subtype;	
}
close INFO;
#print keys(%sampletype);
#my @keys1 = values %sampletype;
#print "$keys1[0]\n";


open OUT,">$out" or die $!;

my %bgm_mutations_p;
my $bgm_count_p = 0;
open BGM,"$bgm_p" or die $!;
while (<BGM>){
	next if /^chrom/;
	my @tmp = split("\t",$_);
	$bgm_count_p += 1;
	# insertion
	my $pos = $tmp[1];
	my ($ref,$alt);
	if (length($tmp[2]) == 1 and length($tmp[3]) > 1){
		$ref = "-";
		$alt = uc(substr($tmp[3],1));
	}elsif (length($tmp[2]) > 1 and length($tmp[3]) == 1){
		$ref = uc(substr($tmp[2],1));
		$alt = "-";
		$pos += 1;
	}else{
		$ref = uc($tmp[2]);
		$alt = uc($tmp[3]);
	}
	my @aa =($tmp[0],$pos,$ref,$alt);
	my $key = join("\t",@aa);
	if (not defined $record{$key}){
		next;
	}
	#print "#$key\n";
	my @bb = @tmp[4..$#tmp];
	$bgm_mutations_p{$key} = \@bb;
}
close BGM;

my %bgm_mutations_t;
my $bgm_count_t = 0;
open BGM1,"$bgm_t" or die $!;
while (<BGM1>){
	next if /^chrom/;
	my @tmp = split("\t",$_);
	$bgm_count_t += 1;
	# insertion
	my $pos = $tmp[1];
	my ($ref,$alt);
	if (length($tmp[2]) == 1 and length($tmp[3]) > 1){
		$ref = "-";
		$alt = uc(substr($tmp[3],1));
	}elsif (length($tmp[2]) > 1 and length($tmp[3]) == 1){
		$ref = uc(substr($tmp[2],1));
		$alt = "-";
		$pos += 1;
	}else{
		$ref = uc($tmp[2]);
		$alt = uc($tmp[3]);
	}
	my @aa =($tmp[0],$pos,$ref,$alt);
	my $key = join("\t",@aa);
	if (not defined $record{$key}){
		next;
	}
	#print "#$key\n";
	my @bb = @tmp[4..$#tmp];
	$bgm_mutations_t{$key} = \@bb;
}

close BGM1;

#my $R = Statistics::R->new();
#$R->start();
#$R->set('z', $pvalue/$bgm_count);
#my $zscore = $R->get('qnorm(1 - (z))');
#my $ad_pvalue = $pvalue/$bgm_count; #BC


for (@lines){
	chomp;
	if (/^chrom/){
		print OUT "$_\ttimes\tpvalue\tadj_pvalue\tmethod\ttag\n";
		next;
	}
	next if /strand_filtered$/;
        my @tmp = split("\t",$_);
        if (scalar @tmp >17){
                my @a = @tmp[0..16];
                @tmp = @a;
        }
        my ($sample,$chrom,$pos,$ref,$alt) = ($tmp[0],$tmp[1],$tmp[2],$tmp[4],$tmp[5]);
        my $key = join("\t",($chrom,$pos,$ref,$alt));
	my %bgm_mutations_cur;
	my $bgm_count_cur;
	print "$sample\n";
	if ($sampletype{$sample} eq "P\n"){
	%bgm_mutations_cur = %bgm_mutations_p;
	$bgm_count_cur = $bgm_count_p;	
	print "Plasma\n";
	}else{
	%bgm_mutations_cur = %bgm_mutations_t;
	$bgm_count_cur = $bgm_count_t;
	print "Tissue\n";
	}
	print "$bgm_count_cur\n";
        if (exists $bgm_mutations_cur{$key}){
		my ($total_count,$alt_count);
		$total_count = $tmp[-4];
		$alt_count = $tmp[-3] - $tmp[-2];
		my ($times,$pvaluet,$aj_pvalue,$method) = filterBGM($key,$total_count,$alt_count,\%bgm_mutations_cur,$bgm_count_cur);
		#if ($times >= 5 && $pvaluet < $pvalue){
		if ($pvaluet < $pvalue){
			push @tmp,($times,$pvaluet,$aj_pvalue,$method,"bgm_pass");
		}else{
			push @tmp,($times,$pvaluet,$aj_pvalue,$method,"bgm_filter");	
		}
	}else{
		push @tmp,("NA","NA","NA","NA","bgm_norecord");
	}
	print OUT join("\t",@tmp),"\n";
}

sub filterBGM{
	my ($key,$total_count,$alt_count,$bgm_mutations,$bgm_count) = @_;
	my %bgm_mutations = %$bgm_mutations;
	my $vaf = $alt_count/$total_count;
	my ($weibull_times,$weibull_pvalue,$z_test_times,$z_test_pvalue,$binom_time,$binom_pvalue) = ("NA","NA","NA","NA","NA","NA");
	my ($aj_weibull_pvalue,$aj_z_test_pvalue,$aj_binom_pvalue) = ("NA","NA","NA");
	my @infos = @{$bgm_mutations{$key}};
	#353903  5       1.41e-05        5/26    1.24115384615385e-05    2.62606523441148e-05    7.463097        0.006864796     0.9176783       0.02800068
	my ($mut,$t) = split("/",$infos[3]);
	if ($mut == 1){
		($binom_time,$binom_pvalue,$aj_binom_pvalue) = binom_test($total_count,$alt_count,$infos[2],$bgm_count);
		return ($binom_time,$binom_pvalue,$aj_binom_pvalue,"binom_test");
	}
	if ($mut <5 || $infos[6] eq "NA" || $infos[7] eq "NA"){
		($z_test_times,$z_test_pvalue,$aj_z_test_pvalue) = z_test($vaf,$infos[4],$infos[5],$bgm_count);
		return ($z_test_times,$z_test_pvalue,$aj_z_test_pvalue,"z_test");
	}else{
		($weibull_times,$weibull_pvalue,$aj_weibull_pvalue) = weibull_test($mut,$t,$vaf,$infos[6],$infos[7],$bgm_count);
		return ($weibull_times,$weibull_pvalue,$aj_weibull_pvalue,"weibull_test");
	}
	
}

sub weibull_test{
	my ($mut,$t,$vaf,$shape,$scale,$bgm_count) = @_;
	my $frac = $mut/$t;
	my $R = Statistics::R->new();
	$R->start();
	$R->set('x', 100 * $vaf);
	$R->set('shape', $shape);
	$R->set('scale', $scale);
	my $WP = 1 - ((1-$frac) + ($frac * $R->get('pweibull(x,shape,scale,TRUE)')));
	my $aj_WP = $WP * $bgm_count;
	$aj_WP = 1 if $aj_WP >1;
	my $mean = $R->get('mean(rweibull(1000,shape=shape,scale=scale))');
	my $times = 100 * $vaf/$mean;
	return ($times,$WP,$aj_WP);
}

sub z_test{
	my ($AF,$meanAF,$stdAF,$bgm_count) = @_;
	#my $AFZ = ($AF - $meanAF) / $stdAF;
	my $R = Statistics::R->new();
	$R->start();
	$R->set('af', $AF);
	$R->set('meanAF',$meanAF);
	$R->set('stdAF',"$stdAF"); # my $zscore = $R->get('qnorm(1 - (z))');
	my $z_pvalue = $R->get('1 - pnorm(af,mean=meanAF,sd=stdAF)');
	my $aj_z_pvalue = $bgm_count * $z_pvalue;
	$aj_z_pvalue =1 if $aj_z_pvalue >1;
	my $z_times = $AF/$meanAF;
	$R->stop();
	return ($z_times,$z_pvalue,$aj_z_pvalue);
}	

sub binom_test{
	my ($total_count,$alt_count,$bgm_vaf,$bgm_count) = @_;
	my $R = Statistics::R->new();
	$R->start();
	$R->set('p', $bgm_vaf);
	$R->set('x',$alt_count);
	$R->set('n',$total_count);
	my $b_pvalue = $R->get('binom.test(x, n, p = p,alternative="greater")$p.value');
	my $aj_b_pvalue = $bgm_count * $b_pvalue;
	$aj_b_pvalue = 1 if $aj_b_pvalue >1;
	$R->stop();
	my $times = ($alt_count/$total_count)/$bgm_vaf;
	return ($times,$b_pvalue,$aj_b_pvalue);
}

