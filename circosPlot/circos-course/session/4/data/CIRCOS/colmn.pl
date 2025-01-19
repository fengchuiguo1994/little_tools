#!/bin/perl

#extract col 1 and 5 from *.ptt files
$i=0;
while(<>)
  {
s#\.\.# #go;
@tab=split(/\s+/, $_);
if ( $i > 2 ) {print "$tab[6] $tab[0] $tab[1]\n"; }
$i++;
}
