dnds.pl
=======

I wrote this script in 2005 to calculate dS/dN ratios from first principles.

It takes as input a file containined pair of aligned sequences in FASTA format.

This distro comes with the standard codon table EGC.0 provided (in 2005) in EMBOSS. If you need a different one, I'm sure you can find it [here](ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz) somewhere.

Put all the files in one directory and run it from there.

```
 perl dnds.pl [options] align.fas
 
 align.fas should contain a pair of codon-aligned sequences
 in fasta format

 dnds strips gaps automatically

 options:
 -t[nn] : desired codon table (emboss EGC.nn)
 -M[c]  : calculation method-
          N[ei] | L[i-Wu-Luo] | P[amilo]
 -R[f]  : transition/tranversion ratio (for -MN)
 ```

MAJ
