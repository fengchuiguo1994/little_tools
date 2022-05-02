#!/usr/bin/perl
#data 2017.03.04
#author xyhuang
#function transform vcf to gwas format
#-------------------------------------------------------------------------
# export package
#-------------------------------------------------------------------------
use strict;
use warnings;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

#---------------------------------------------------------------------------
# get options
#---------------------------------------------------------------------------
my($vcf,$value,$DP);
GetOptions(
  "vcf=s"  => \$vcf,
  "q=i"  => \$value,
  "DP=i"  => \$DP,
  "help|?" => \&USAGE,
) or &USAGE;

&USAGE unless (defined $vcf);
my $tmp = $1 if($vcf=~/(.+?)\.vcf/);
my $info ||= "$tmp.inform";
my $num ||= "$tmp.numeric";
$value ||= 200;
&log_current_time("$Script start...");

#----------------------------------------------------------------------------
# load input file and data processing
#----------------------------------------------------------------------------
open IN,"<$vcf" or die "$vcf is not a vcf file";
open INFO,">$info" or die "$info can't available";
print INFO "SNP\tChromosome\tPosition\n";
open NUM,">$num" or die "$num can't available";
my $number;
while(<IN>){
	next if(/^##/);
	chomp;
	my $line = $_;
	if(/^#/){
		my @tmp=split/\t/,$line;
		$number = $#tmp-8;
		$DP ||= $number*4;
		my @arr;
		$arr[0] = 'taxa';
		push @arr,@tmp[9..$#tmp];
		print NUM (join "\t",@arr)."\n" ;
	}
	else{
		my @tmp=split/\t/,$line;
		next if($tmp[5] < $value);
		if($tmp[7] =~ /;DP=(\d+);/){
			next if($1 < $DP);
		}
		$tmp[0]=~s/[cC][Hh][Rr]//g;
		my $id="$tmp[0]_$tmp[1]";
		my @array;
		push @array,$id;
		my $flag0 = 0;
		my $flag2 = 0;
		foreach(9..$#tmp){
			if($tmp[$_]=~/(.)\/(.)/){
				if($1 ~~ $2){
					if($1 eq '.'){
						push @array,"0";
						$flag0++;
					}
					elsif($1==0){
						push @array,0;
						$flag0++;
					}
					else{
						push @array,2;
						$flag2++;
					}
				}
				else{
					push @array,1;
				}
			}
		}
		next if($flag0/$number < 0.1 || $flag2/$number < 0.1);
		print NUM (join "\t",@array)."\n" ;
		print INFO "$id\t$tmp[0]\t$tmp[1]\n";
	}
}

#----------------------------------------------------------------------------------
# print result and end of work
#----------------------------------------------------------------------------------
&log_current_time("$Script end...");
my $run_time=time()-$BEGIN_TIME;
print "$Script run time :$run_time\.s\n";

#----------------------------------------------------------------------------------
# function
#----------------------------------------------------------------------------------
sub log_current_time {
	my ($info) = @_;
	my $curr_time = &date_time_format(localtime(time()));
	print "[$curr_time] $info\n";
}

sub date_time_format {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE{
	my $usage=<<"__USAGE__";
#-----------------------------------------------------------
Program:$Script
Version:$version
Contact:1182768992\@qq.com
Data:2017-03-04
Function:transform vcf to gwas file
USAGE:
	--vcf	<STR>   input vcf file,output loci file and genetype file [Must]
	--q		<INT>   filter by qvalue	[Default:200]
	--DP	<INT>   all sample depth	[Default:sample number * 4]
	--help          show the docment and exit
Example:
perl $Script --vcf merge.vcf
#---------------------------------------------------------
__USAGE__
   print $usage;
   exit;
}
