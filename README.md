# Radon

Module to perform Direct and Inverse Discrete Radon Transforms in PERL
using PDL. Works on square images with size 2^n.

Based on *Discrete Radon transform has an exact, fast inverse and
generalizes to operations other than sums along lines*, William H. Press
PNAS (2006) 103 (51) 19249-19254
https://doi.org/10.1073/pnas.0609228103.

For its use see the example at testradon.pl.

The plots in 'reconstructed.pdf' were obtained through the command
`$ ./testradon.pl --image=test-pgm --iterations=5`

## Installation

This is not yet a well formed package. Just unpack and copy the files
to any convenient place. If the file *Radon.pm* is not in your @INC
array, add the appropriate `use lib` to your program before you `use
Radon`. To test the example you might have to adjust that line in
*testradon.pl* and provide the path to the test image.

## Author

W. Luis Mochán mochan@fis.unam.mx

## Acknowledgment

This work was partially supported by DGAPA-UNAM under grant IN111119.

## Licence

This software is copyright (c) 2020 by W. Luis Mochan.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.
