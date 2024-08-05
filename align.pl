use strict;
use warnings;
use threads;
use threads::shared;
use Cwd qw (abs_path);
use FindBin qw ($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;
use Getopt::Long;
use Data::Dumper;
use Littletools;

my $BEGIN_TIME=time();
my $version="1.0.0";

my($data,$info,$outdir);
GetOptions(
  "data=s"  => \$data,
  "info=s"  => \$info,
  "outdir=s"  => \$outdir,
  "help|?" => \&USAGE,
 ) or &USAGE;

&USAGE unless (defined $data and defined $info);
sub USAGE{
	my $usage = <<"__USAGE__";
#--------------------------------------------------------------
  Program:$Script
  Version:$version
  Contact:1182768992\@qq.com
     Date:2017-2-27
 Function:a pipeline for DNA data analysis
    USAGE:
		--data		<STR>	data input file[Must]
		--info		<STR>	info input file[Must]
		--outdir	<STR>	output dictionary[Optional default:./]
		--log		<STR>	output the log file[Optional default:./run.log]
		--help/h	<STR>	show the docment and exit
  Example:
		perl $Script --data data.txt --info info.txt --outdir ./
#--------------------------------------------------------------
__USAGE__
	print $usage;
	exit;
}

print &log_current_time("$Script start ...")."\n";
$outdir ||= "./";
$outdir = abs_path($outdir);
system "mkdir -p $outdir" unless (-d $outdir);
print &log_current_time("loading data and info start ...")."\n";
my $start = time();
my %data = %{&readfile2($data)};
my @sample;
my @data = keys %data;
foreach(@data){
	if(scalar @{$data{$_}} > 2 || scalar @{$data{$_}} < 1){
		print "please check the data.info ,the RNA-Seq data is single_end or pair_end,why the data.info has over two data or less than one\n";
		exit;
	}
	push @sample,$_ foreach(@{$data{$_}});
}
my %info = %{&readfile1($info)};
my $use = time()-$start;
print &log_current_time("loading data and info finished [run time :$use\.s]")."\n";

print &log_current_time("start the pipeline ...")."\n";
my (%thr,%length);
$info{threads} ||= 1;

print &log_current_time("this step is building genome index ...")."\n";
if($info{genomeidx} ne 'none'){
	print "there are genome index file , it will be used\n";
}
else{
	print "there are not genome index file , using genome file to build index\n";
	my $idxout = "$outdir/genomeidex";
	system "mkdir -p $idxout" unless (-d $idxout);
	print &log_current_time("building genome index start...")."\n";
	my $start = time();
	system("$info{star} --runMode genomeGenerate --genomeDir $idxout --genomeFastaFiles $info{genome} --sjdbGTFfile $info{annotation} --sjdbOverhang $info{readslength} --runThreadN $info{threads}");
	my $use = time()-$start;
	print &log_current_time("building genome index finish: $info{star} --runMode genomeGenerate --genomeDir $idxout --genomeFastaFiles $info{genome} --sjdbGTFfile $info{annotation} --sjdbOverhang $info{readslength} --runThreadN $info{threads} [run time :$use\.s]")."\n";
	$info{genomeidx} = $idxout;
}

print &log_current_time("this step is align ...")."\n";
my $align = "$outdir/align";
system "mkdir -p $align" unless (-d $align);
foreach(@data)
{
	my @tmp = @{$data{$_}};
	my $start = time();
	print &log_current_time("STAR align $_ start...")."\n";
	system ("STAR --runThreadN $info{threads} --genomeDir $info{genomeidx} --readFilesIn $tmp[0] $tmp[1] --readFilesCommand zcat --sjdbGTFfile $info{annotation} --outFilterMultimapNmax 12 --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFileNamePrefix $align/$_ --outSAMattrIHstart 0 --outSAMmapqUnique 60 --outSAMmultNmax 10 --outSAMattributes NH HI AS nM NM MD jM jI XS --outSAMtype BAM Unsorted --chimSegmentMin 35 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic");
	my $use = time()-$start;
	print &log_current_time("STAR align $_ finish: STAR --runThreadN $info{threads} --genomeDir $info{genomeidx} --readFilesIn $tmp[0] $tmp[1] --readFilesCommand zcat --sjdbGTFfile $info{annotation} --outFilterMultimapNmax 12 --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFileNamePrefix $align/$_ --outSAMattrIHstart 0 --outSAMmapqUnique 60 --outSAMmultNmax 10 --outSAMattributes NH HI AS nM NM MD jM jI XS --outSAMtype BAM Unsorted --chimSegmentMin 35 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic [run time :$use\.s]")."\n";

}

print &log_current_time("$Script end¡­¡­")."\n";
my $run_time=time()-$BEGIN_TIME;
print "$Script run time :$run_time\.s\n";
