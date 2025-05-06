use strict;
use warnings;

my ($flag,%hashtmp);
while(<>){
	chomp;
	my @tmp = split /\t/,$_;
#	my $spe = (split /\./,$tmp[16])[0];
	my $spe = $tmp[48];
	if(defined $flag && $flag ne $tmp[0]){
		print "$flag\t".(sort {$hashtmp{$b} <=> $hashtmp{$a}} keys %hashtmp)[0]."\n";
		undef %hashtmp;
	}
	$flag = $tmp[0];
	$hashtmp{$spe} += 1;
}
print "$flag\t".(sort {$hashtmp{$b} <=> $hashtmp{$a}} keys %hashtmp)[0]."\n";
