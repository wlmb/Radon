#!/usr/bin/env perl
# Perform and display a direct and inverse discrete Radon transform
#
use warnings;
use strict;
use lib "."; # Adjust as needed to pick up Radon.pm if uninstalled
use feature qw(say);
use IO::Prompter;
use PDL;
use PDL::NiceSlice;
use PDL::IO::Pic;
use PDL::Image2D;
use PDL::Graphics::Gnuplot;
use Getopt::Long;
use Radon;

my $name;
my $iterations=0;
GetOptions("image=s"=>\$name, "iterations=i"=>\$iterations)
    or usage('Bad options');
usage("Missing image filename") unless defined $name;
usage("Can't read $name") unless -r $name;
my $im=rpic($name);
usage("Image should be monochromatic") unless $im->ndims==2;
my $N=$im->dim(0);
usage("Image should be square") unless $im->dim(1)==$N;
usage("Size should be power of 2") if $N&($N-1);
my ($N0,$N1, $N2, $N3, $N4)=map {$_*$N} (0..4);
my ($a, $b, $c, $d)=Radon::radonD($im);
my $rall=Radon::radonA($a, $b, $c, $d);
my $im0=Radon::radonI($a, $b, $c, $d, 0);
my $imN=Radon::radonI($a, $b, $c, $d, $iterations);
my $w1=gpwin('qt', size=>[12,9]);
$w1->multiplot(Layout=>[2,2]);
$w1->plot({title=>'Original'}, with=>'image',$im);
$w1->plot({title=>'Transform'}, with=>'image',$rall);
$w1->plot({title=>'Reconstructed: 0 iterations'}, with=>'image',$im0);
$w1->plot({title=>"Reconstructed: $iterations iterations"}, with=>'image',$imN);
$w1->multi_end;

sub usage {
    say $_[0];
    say <<'FIN';
Usage:
   ./testradon.pl --image=path --iterations=n

Tests the package Radon by calculating the discrete radon transform of
the image at path and the inverse transform using n refinement
iterations. The image should be monochrome, square and its size should
be a power of 2.
FIN
    exit 1;
}
