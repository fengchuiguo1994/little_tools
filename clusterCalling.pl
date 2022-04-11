use strict;
use Data::Dumper;

# chr1  start1  end1    chr2    start2  end2    <interactionCount> <... any other info ...>
# sort -k1,1 -k4,4 -k2,2n -k3,3n -k5,5n -k6,6n | perl clusterCalling.pl > result
# or 
# sort -k1,1 -k4,4 -k2,2n -k3,3n -k5,5n -k6,6n tmp.file 
# perl clusterCalling.pl tmp.file > result

my $chrom = "";
my $end;
my @part;

while(<>){
    chomp;
    if ($_ =~ /^\s*$/) { # continue blank records
        next;
    }
    my @tmp = split(/\s+/,$_);
    if (scalar(@tmp) == 6 ) { # default interaction count
        push(@tmp,1);
    }

    push(@tmp,0); # mark

    if ($chrom ne "" && ($chrom ne $tmp[0] or $end < $tmp[1])){ # has no overlap
        cluster(@part); # cluster calling
        undef @part;
    }
    $chrom = $tmp[0];
    push(@part,\@tmp);
    $end = $tmp[2];
}
cluster(@part); # cluster calling

sub cluster{
    my @rawInteraction;
    my @newInteraction = @_;
    while (scalar(@rawInteraction) != scalar(@newInteraction)) {
        @rawInteraction = @newInteraction;
        foreach my $i (0..$#rawInteraction-1){
            if ($rawInteraction[$i][[-1]] == 1) { # has been merged
                next;
            }
            foreach my $j ($i+1..$#rawInteraction){
                if ($rawInteraction[$j][[-1]] == 1) { # has been merged
                    next;
                }
                if (bedpeOverlap($rawInteraction[$i],$rawInteraction[$j])){
                    mergeBedpe($rawInteraction[$i],$rawInteraction[$j]);
                }
            }
        }
        undef @newInteraction;
        foreach my $i (0..$#rawInteraction){
            if ($rawInteraction[$i][-1] == 0){
                push(@newInteraction,$rawInteraction[$i]);
            }
        }
    }
    foreach my $i (0..$#newInteraction){
        my @tmp = @{$newInteraction[$i]};
        pop(@tmp);
        print join("\t",@tmp)."\n";
    }
}


sub bedpeOverlap{
    my $tmp1 = shift;
    my $tmp2 = shift;
    my @bed1 = @{$tmp1};
    my @bed2 = @{$tmp2};
    if ($bed1[0] eq $bed2[0] && $bed1[3] eq $bed2[3]){
        my $tmp = bedOverlap($bed1[1],$bed1[2],$bed2[1],$bed2[2]) && bedOverlap($bed1[4],$bed1[5],$bed2[4],$bed2[5]);
        return bedOverlap($bed1[1],$bed1[2],$bed2[1],$bed2[2]) && bedOverlap($bed1[4],$bed1[5],$bed2[4],$bed2[5]);
    }
    return 0;
}

sub bedOverlap{
    my $s1 = shift;
    my $e1 = shift;
    my $s2 = shift;
    my $e2 = shift;
    return $e2 > $s1 && $s2 < $e1; # or return $e2 >= $s1 && $s2 <= $e1;
}

sub mergeBedpe{
    my $bed1 = shift;
    my $bed2 = shift;    
    if (${$bed1}[1] > ${$bed2}[1]){
        ${$bed1}[1] = ${$bed2}[1]
    }
    if (${$bed1}[2] < ${$bed2}[2]){
        ${$bed1}[2] = ${$bed2}[2]
    }
    if (${$bed1}[4] > ${$bed2}[4]){
        ${$bed1}[4] = ${$bed2}[4]
    }
    if (${$bed1}[5] < ${$bed2}[5]){
        ${$bed1}[5] = ${$bed2}[5]
    }
    ${$bed2}[-1] = 1;
    ${$bed1}[6] += ${$bed2}[6]; 
}


