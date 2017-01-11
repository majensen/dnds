dnds.pl
=======

I wrote this script in 2005 to calculate dS/dN ratios from first principles.

It takes as input a file containing pair of aligned sequences in FASTA format.
(The fasta can contain multiple such pairs: seq 1+2, 3+4, etc. will be analyzed.)

This distro comes with the standard codon table EGC.0 provided (in 2005) in EMBOSS. If you need a different one, I'm sure you can find it [here](ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz) somewhere.

Put all the files in one directory and run it from there.

Do the following:
```
perl dnds.pl -t0 -MN test.fas
```
and you should get back
```
seq1    seq2    ds      dnv     dn      dsv
>MG1655_m56_94; >SenLT2_v1_90773;       0.671179        3.87855e-05     0.0759739       0.00221391
```

You can choose between Nei-Gojobori (1986) and Li-Wu-Luo (1985) methods.


```
 perl dnds.pl [options] align.fas
 
 align.fas should contain a pair of codon-aligned sequences
 in fasta format

 dnds strips gaps automatically

 options:
 -t[nn] : desired codon table (emboss EGC.nn)
 -M[c]  : calculation method-
          N[ei] | L[i-Wu-Luo]
 -R[f]  : transition/tranversion ratio (for -MN)
 ```

MAJ
