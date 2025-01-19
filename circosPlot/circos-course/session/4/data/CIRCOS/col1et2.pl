#!/bin/perl

#extract col 1 and 2 from a file

while(<>)
{
    @tab=split(/\s+/, $_);

    print"$tab[0]\t$tab[1]\n";
}
