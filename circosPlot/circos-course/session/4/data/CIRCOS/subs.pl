#!/bin/perl

#substitute gene by chrname

#SACE_COORD, CAGL_COORD
while(<>)
  {
@tab=split(/\s+/, $_, 2);
if ( m/^YA/ ) { print "$tab[0]\tsace-a\t$tab[1]";}
if ( m/^YB/ ) { print "$tab[0]\tsace-b\t$tab[1]";}
if ( m/^YC/ ) { print "$tab[0]\tsace-c\t$tab[1]";}
if ( m/^YD/ ) { print "$tab[0]\tsace-d\t$tab[1]";}
if ( m/^YE/ ) { print "$tab[0]\tsace-e\t$tab[1]";}
if ( m/^YF/ ) { print "$tab[0]\tsace-f\t$tab[1]";}
if ( m/^YG/ ) { print "$tab[0]\tsace-g\t$tab[1]";}
if ( m/^YH/ ) { print "$tab[0]\tsace-h\t$tab[1]";}
if ( m/^YI/ ) { print "$tab[0]\tsace-i\t$tab[1]";}
if ( m/^YJ/ ) { print "$tab[0]\tsace-j\t$tab[1]";}
if ( m/^YK/ ) { print "$tab[0]\tsace-k\t$tab[1]";}
if ( m/^YL/ ) { print "$tab[0]\tsace-l\t$tab[1]";}
if ( m/^YM/ ) { print "$tab[0]\tsace-m\t$tab[1]";}
if ( m/^YN/ ) { print "$tab[0]\tsace-n\t$tab[1]";}
if ( m/^YO/ ) { print "$tab[0]\tsace-o\t$tab[1]";}
if ( m/^YP/ ) { print "$tab[0]\tsace-p\t$tab[1]";}

if ( m/^CAGL0A/ ) { print "$tab[0]\tcagl-a\t$tab[1]"; }
if ( m/^CAGL0B/ ) { print "$tab[0]\tcagl-b\t$tab[1]"; }
if ( m/^CAGL0C/ ) { print "$tab[0]\tcagl-c\t$tab[1]"; }
if ( m/^CAGL0D/ ) { print "$tab[0]\tcagl-d\t$tab[1]"; }
if ( m/^CAGL0E/ ) { print "$tab[0]\tcagl-e\t$tab[1]"; }
if ( m/^CAGL0F/ ) { print "$tab[0]\tcagl-f\t$tab[1]"; }
if ( m/^CAGL0G/ ) { print "$tab[0]\tcagl-g\t$tab[1]"; }
if ( m/^CAGL0H/ ) { print "$tab[0]\tcagl-h\t$tab[1]"; }
if ( m/^CAGL0I/ ) { print "$tab[0]\tcagl-i\t$tab[1]"; }
if ( m/^CAGL0J/ ) { print "$tab[0]\tcagl-j\t$tab[1]"; }
if ( m/^CAGL0K/ ) { print "$tab[0]\tcagl-k\t$tab[1]"; }
if ( m/^CAGL0L/ ) { print "$tab[0]\tcagl-l\t$tab[1]"; }
if ( m/^CAGL0M/ ) { print "$tab[0]\tcagl-m\t$tab[1]"; }

if ( m/^ZYRO0A/ ) { print "$tab[0]\tzyro-a\t$tab[1]"; }
if ( m/^ZYRO0B/ ) { print "$tab[0]\tzyro-b\t$tab[1]"; }
if ( m/^ZYRO0C/ ) { print "$tab[0]\tzyro-c\t$tab[1]"; }
if ( m/^ZYRO0D/ ) { print "$tab[0]\tzyro-d\t$tab[1]"; }
if ( m/^ZYRO0E/ ) { print "$tab[0]\tzyro-e\t$tab[1]"; }
if ( m/^ZYRO0F/ ) { print "$tab[0]\tzyro-f\t$tab[1]"; }
if ( m/^ZYRO0G/ ) { print "$tab[0]\tzyro-g\t$tab[1]"; }



  }
