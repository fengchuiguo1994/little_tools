#!/bin/perl

#from intrspecies duplication: remove identical lines.
#SACESACE, CAGLCAGL and ZYROZYRO may include multiple hits between identical
#query and hits sequences

#this script is meant for considering the same order in the original file
#otherwise consider : sort -u CAGLCAGL

$IN=@ARGV[0];  # exp CAGLCAGL
open(IN,"$IN");
$SQ=$SH="";
while(<IN>)
  {

@tab=split(/\s+/, $_);
if ($tab[0] eq $SQ && $tab[1] eq $SH ) { next; }
$SQ=$tab[0];
$SH=$tab[1];

print $_;

  }

close(IN);
