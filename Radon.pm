package Radon;
our $VERSION=0.001;
use strict;
use warnings;

use base 'Exporter';
our @EXPORT_OK=qw(radonD radonI radonA);

use List::Util;
use PDL;
use PDL::NiceSlice;
use PDL::Image2D;


sub radonD { #discreet radon transform
    my $image=shift;
    die "Image should be square" unless $image->dim(0)==$image->dim(1);
    my $nn=$image->dim(0);
    die "Image should have more than one pixel" if $nn<2;
    die "Image size should be power of 2" if $nn & $nn-1;
    my ($a, $b, $c, $d)= map {zeroes($nn,2*$nn)} (1..4);
    $a(:,$nn:-1).=$image;
    $b(:,$nn:-1).=$image->transpose;
    $c(:,$nn:-1).=$image->transpose->(:,-1:0);
    $d(:,$nn:-1).=$image->(-1:0,:);
    return map {$_=xform($_)} ($a, $b, $c, $d); #transform
}

sub xform { #Transform one octant.
    my $a=shift; #Sector. Initial value is image at top, zeroes at bottom
    my $N=$a->dim(0);
    my $N2=$a->dim(1);
    my $downstep=$N; #somehow count down
    my $upstep=1;
    while($downstep!=1){
	#s from 0 to upstep-1 counts angles.
	#n from 0 to downstep-1 counts strips
	#h from 0 to n2-1 counts heights
	#p counts parity (even-odd)
	my $r=$a->copy->reshape($upstep, $downstep, $N2);#s,n,h one plane per s value
	my $Nslices=$downstep>>1;
	my $left=$r->(:,0:-1:2,:); #s,n/2,h left parts
	my $rawright=$r->(:,1:-1:2,:); #s,n/2,h right parts
	my $right=$rawright->range([map {[$_,0,$_]} 0..$upstep-1],
				   [0,$Nslices,$N2], 'p');
	my $even=$left+$right;
	my $odd=$left+$right->mv(-1,0)->rotate(-1)->mv(0,-1);
	#intercalate even and odd results
	my $res=pdl($even, $odd)->mv(-1,0)->reshape($N,$N2); #s,n/2,h,p->p,s,n/2,h->psn,h->s,h
	$upstep<<=1; #multiply by 2
	$downstep>>=1; #divide by 2
	$a=$res;
    }
    return $a;
}

sub backproject { #simple back projection (sum of lines through point)
    my ($a, $b, $c, $d)=@_;
    my $N=$a->dim(0);
    die "Dimension should be square" if $N&$N>>1;
    die "Wrong dims" unless List::Util::all
    {$_->dim(1)==2*$_->dim(0) && $_->dim(0) == $N}
    ($a, $b, $c, $d);
    my $back=backpiece($a) + backpiece($b)->transpose +
	backpiece($c)->(:,-1:0)->transpose + backpiece($d)->(-1:0,:);
    return $back/(4*($N-1));
}

sub backpiece { #backproject one quadrant 0:45, 45_90, -90_-45, or -45:0
    my $a=shift; #Sector. Initial value is image at top, zeroes at bottom
    my $N=$a->dim(0);
    my $N2=$a->dim(1);
    my $downstep=$N; #somehow count down
    my $upstep=1;
    while($downstep!=1){
	#s from 0 to downstep-1 counts angles.
	#n from 0 to upstep-1 counts images
	#h from 0 to n2-1 counts heights
	#p counts parity (even-odd)
	#s' from 0 to downstep/2-1 counts new angles
	my $Ns1=$downstep>>1; #number of s'
	my $r=$a->copy->reshape($downstep, $upstep, $N2); #s,n,h
	my $Nslices=$upstep<<1;
	my $even=$r->(0:-1:2,:,:); #s',n,h; s' is floor(s/2)
	my $odd =$r->(1:-1:2,:,:); #s',n,h
	my $left=$even+$odd; #s',n,h
	my $rawright=$even+$odd->mv(-1,0)->rotate(1)->mv(0,-1); #s',n,h
	my $right=$rawright->range([map {[$_,0,-$_]} 0..$Ns1-1],
				   [0,$upstep,$N2], 'p');
	#$right*=($right->xvals<$right->yvals+1);
	my $res=pdl($left, $right)->mv(-1,1)->reshape($N, $N2);	#s',n,h,lr->s',lr,n,h->s'n',lr
	$downstep>>=1;
	$upstep<<=1;
	$a=$res;
    }
    return $a->(:,$N:-1);
}

sub radonI { #inverse of discreet radon transform
    my ($a, $b, $c, $d, $iterations)=@_;
    my $N=$a->dim(0);
    die "Dimension should be square" if $N&$N>>1;
    die "Wrong dims" unless List::Util::all
    {$_->dim(1)==2*$_->dim(0) && $_->dim(0) == $N}
    ($a, $b, $c, $d);
    my $f=recursive_inverse($a, $b, $c, $d);
    for(1..$iterations){
	my ($a1, $b1, $c1, $d1)=radonD($f);
	my ($resa, $resb, $resc, $resd)=($a-$a1, $b-$b1, $c-$c1, $d-$d1); #residuals
	my $rf=recursive_inverse($resa, $resb, $resc, $resd); #residual image
	$f=$f+$rf; #correct
    }
    return $f;
}

sub recursive_inverse { #recursively invert radon transformation through all scales
    my ($a, $b, $c, $d)=@_;
    my $N=$a->dim(0);
    return $a->(0,1) if $N==1;
    my ($sa, $sb, $sc, $sd)= #restrict
	map {($_->(0:-1:2,0:-1:2)+$_->(0:-1:2,1:-1:2))/4} ($a, $b, $c, $d);
    my $f2=recursive_inverse($sa, $sb, $sc, $sd); #recursion at half resolution
    my $lowf=$f2->(*2,:,*2,:)->reshape($N,$N); #low resolution
    my ($lowa, $lowb, $lowc, $lowd)=radonD($lowf); # radon of low resolution
    my ($resa, $resb, $resc, $resd)=($a-$lowa, $b-$lowb, $c-$lowc, $d-$lowd); #residuals
    my $res=backproject($resa, $resb, $resc, $resd); #residual image
    my $kern=pdl([[-1/16, -1/8, -1/16],
		  [-1/8, 3/4, -1/8],
		  [-1/16, -1/8, -1/16]]); #high frequency kernel
    my $reshp=$res->conv2d($kern); #high pass filter
    my $f=$lowf+$reshp;
    return $f;
}

sub radonA { #assemble the four pieces in one for display
    #Note that height increases with second index but has different
    #meaning in different quadrants. Rise goes towards the right,
    #left, right left according to quadrant.
    #
    my ($a, $b, $c, $d)=@_;
    my $N=$a->dim(0);
    my ($N0,$N1, $N2, $N3, $N4)=map {$_*$N} (0..4);
    my $res=zeroes($N4,$N2);
    $res(0:$N1-1,:).=$c(:,-1:0);
    $res($N1:$N2-1,:).=$d(-1:0,:);
    $res($N2:$N3-1,:).=$a;
    $res($N3:$N4-1,:).=$b(-1:0,-1:0);
    return $res;
}
