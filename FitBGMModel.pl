use strict;
use warnings;
use Statistics::Descriptive;
use Statistics::R;
use Getopt::Long;

my ($in,$out,$filter_sample,$filter_max,$help);

GetOptions(
  "i|in=s"=>\$in,
  "o|out=s"=>\$out,
  "f|filter_sample=s"=>\$filter_sample,
  "m|filter_max"=>\$filter_max,
  "h|help"=>\$help,
);

usage() && exit if $help|| !$in||!$out;

sub usage{
	print"Usage: perl $0 -i <infile> -o <outfile> [-f filter_sample]\n";
}

open IN,"$in" or die $!;
my $header = <IN>;
chomp $header;
my @h = split("\t",$header);
my $filtered_index = "NA";
my $total_sample = $#h-7+1;
for (7..$#h){
	if ( $filter_sample && $h[$_] eq $filter_sample){
		$filtered_index = $_;
		$total_sample -= 1;
	}
}

my $R = Statistics::R->new();
$R->startR;
#$R->run(q`.libPaths(c("/share/work3/capsmart/pipeline/capSMART/CAPcSMART/Rlib","/share/work3/wangrr/local/Rlib","/share/work2/liuyan/RlocalLib","/home/capsmart/R/x86_64-pc-linux-gnu-library/3.4","/share/work3/wangrr/local/anaconda2/lib/R/library"))`);
$R->run(q`library(fitdistrplus)`);

open OUT,">$out" or die $!;
print OUT "chrom\tpos\tref\talt\ttotal_count\talt_count\tvaf\tMut/Nonmut\tMeanAF\tStdAF\tW_Shape\tW_Scale\tW_Corr\tW_Pval\n";

while (<IN>){
	chomp;
	next if /^chrom/;
	my @tmp = split("\t",$_);
	next if $tmp[5] <2;
	next if $tmp[3] eq "N";
	my @olist = @tmp[0..6];
	my @vafs;
	for my $i (7..$#tmp){
		if ($filtered_index && $i eq $filtered_index){
			#print "$i\n";
			next;
		}else{
			my $aa = (split(",",$tmp[$i]))[-1];
			push @vafs,$aa;
		}
	}
	my $count = 0;
	for (@vafs){
		if ($_ > 0){
			$count += 1;
		}
	}
	if ($count <2){
		my $mut_vs_nomut = "$count/$total_sample";
		push @olist,($mut_vs_nomut,$tmp[6],"NA","NA","NA","NA","NA");
		print OUT join("\t",@olist),"\n";
		next;
	}
	my ($meanAF,$stdAF) = gauss_fit(@vafs);
	my ($W_Shape,$W_Scale,$W_Corr,$W_Pval) = weibull_fit(@vafs);
	my $mut_vs_nomut = "$count/$total_sample";
	push @olist,($mut_vs_nomut,$meanAF,$stdAF,$W_Shape,$W_Scale,$W_Corr,$W_Pval);
	print OUT join("\t",@olist),"\n";
}


sub gauss_fit{
	my @vs = @_;
	my $vafobj = Statistics::Descriptive::Full->new();
	$vafobj->add_data(@vs);
	return ($vafobj->mean(),$vafobj->standard_deviation());
}

sub weibull_fit{
	my @temp = @_;
	my $shape = "NA";
	my $scale = "NA";
	my $cor = "NA";
	my $cpval = "NA";
	$R->set('x', \@temp);
	$R->run(q`x <- x[which(x>0)]`);
	my $len = $R->get('length(unique(x))');
	if($len>2){
		$R->run(q`w <- NA`);
		$R->run(q`tryCatch({w<-suppressWarnings(summary(fitdist(x*100, "weibull"))$estimate)}, warning = function(war){ b<-1 }, error = function(err){ b<-1 }, finally = { print(w) })`);
		if($R->get('w[1]') ne "NA"){ #if parameter estimation successful
			$shape = $R->get('w[1]');
			$scale = $R->get('w[2]');
        	    	$R->run(q`w <- qqplot(x*100, rweibull(100, w[1], w[2]))`);
	        	$R->run(q`co <- cor.test(w$x, w$y)`);
	        	$cor = $R->get('co$estimate');
		        $cpval = $R->get('co$p.value');
		 }
	}
	return ($shape,$scale,$cor,$cpval);
}


