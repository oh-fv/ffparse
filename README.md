# ffparse.py

Parses the output of a Firefly[1] XMCQDPT2 calculation.

* prints information about the active space
* prints MP2 states in energetic order with their contribution from the CI states
* prints the CAS excitations and their contribution to the CI states
* if the calculation was set up using symmetry, the symmetry of the orbitals and excitations is included in the output

The code is very ugly and not guaranteed to work with your output file. Please adapt it to your needs.

## Usage
`$ ./ffparse.py firefly.out > firefly-results.txt`

## Examples
Example input and output files are included in the `examples` folder.

## Dependencies
* O. Hüter, [molsym](https://github.com/oh-fv/molsym) package, Analytic point group algebra for molecular symmetry operations.

## References
1. Alex A. Granovsky, Firefly version 8, [http://classic.chem.msu.su/gran/gamess/index.html](http://classic.chem.msu.su/gran/gamess/index.html)