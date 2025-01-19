#!/bin/perl

#create link file 
#3SPEC_newCOORD  SACESACE (all thits)

$C=@ARGV[0]; # 3SPEC_newCOORD file
$H=@ARGV[1]; #file including hits *SACESACE, SACECAGL,...)
open(IN,"$C");
while(<IN>)
  {

chomp($_);
@tab=split(/\s+/, $_, 2);
$COORD{$tab[0]} = $tab[1];

  }
close(IN);

open(H,"$H");
while(<H>)
  {

chomp($_);
@tab=split(/\s+/, $_ );

print "$COORD{$tab[0]} $COORD{$tab[1]}\n";

  }
