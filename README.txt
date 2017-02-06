

===========================================================
SweeD
===========================================================


Copyright (C) 2012 Pavlos Pavlidis and Nikolaos Alachiotis 


SweeD implements a maximum-likelihood approach to detect selective sweeps 
using the site frequency spectrum.

The sequential version of SweeD represents an algorithmically improved re-design
of the SweepFinder tool by Nielsen et al. (2005). 

For more details see:
R. Nielsen, S. Williamson, Y. Kim, M.J. Hubisz, A.G. Clark, and C. Bustamante, 
"Genomic scans for selective sweeps using SNP data", Genome Research, pages 1566-1575, 2005.


Linux Platforms
---------------

Compile:
	make -f Makefile.gcc

Execute:
	./SweeD -name TEST -input TEST.SF -grid 10000


Windows Platforms
-----------------

No Windows version is yet available.


Change Log
----------
March 2014: SweeD v3.2.4 (it can read the ms-like output of msABC)

October 2013: SweeD v3.2.1 (fixed bug in VCF format associated with ALT state as '.')

January 2013:	SweeD v3.1 (fixed bug in VCF parser associated with the handling of missing data)

November 2012:	SweeD v3.0 (analytical calculation of SFS)

November 2012:	SweeD v2.6 (fixed bug in multiple SFS input file)

October 2012:	SweeD v2.5 (checkpointable)

October 2012:	SweeD v2.0 (scalable for more than 1000 sequences, parallel)

September 2012:	SweeD v1.0


Contact
-------

Pavlos Pavlidis - pavlidisp@gmail.com
Nikos Alachiotis - n.alachiotis@gmail.com


This program is free software; you may redistribute it and/or modify its
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option)
any later version.
 
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.
