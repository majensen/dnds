use Getopt::Std;
use Cdn;
# note: convert U to T, lower to upper case
# dnds: calculate dn, ds for two codon-aligned sequences
#       using various methods

# perl dnds [options] align.fas
# 
# align.fas should contain a pair of codon-aligned sequences
# in fasta format
#
# dnds strips gaps automatically
#
# options:
# -t[nn] : desired codon table (emboss EGC.nn)
# -M[c]  : calculation method-
#          N[ei] | L[i-Wu-Luo]
# -R[f]  : transition/tranversion ratio (for -MN)


getopts("t:M:R:");

#$emb_loc = "c:/cygwin/usr/local/share/EMBOSS/data";
$emb_loc = ".";

#defaults and checking
$opt_t = 0 unless ($opt_t);
$opt_M = "N" unless ($opt_M);
die "dnds: -M method '$opt_M' not recognized" if ($opt_M !~ /^[LNP]/);
if ($opt_M =~ /^N/) {
    $opt_R = 0.5 unless $opt_R;
}
else {
    warn "dnds: -R ignored for -M$opt_M\n" if $opt_R;
}

die "dnds: need alignment file" unless $ARGV[0];
open ALN, $ARGV[0] or die "dnds: alignment file problem:$!";

# build codon and degenerate sites table

die "dnds: can't find codon table" unless (-e "$emb_loc/EGC.$opt_t");

# build codon tables
$cdntab = read_codon_table("$emb_loc/EGC.$opt_t");
$cdnlist = pack "(a3)*", keys %{$cdntab};
$sites = count_sites( $cdnlist, $cdntab );
@sitab{unpack ("(a3)*",$cdnlist)} = unpack("(a3)*",$sites);
$sitab = \%sitab;  # for consistency below
$pthtab = path_table( $cdntab, $sitab );
$cdntab->{"---"} = "---";
$sitab{"---"} = "---";

# get sequences
$i=-1;
while (<ALN>) { 
    chomp;
    if (/>/) {
	$i++;
	$nm[$i] = $_;
	$seq[$i] = "";
    } 
    else {
	$seq[$i] .= uc;
    }
}
close ALN;

die "dnds: only one sequence in $ARGV[0]" unless $i >= 1;
$i-- unless $i % 2; 


for $j (0..$i/2) {
    warn ("dnds: sequences ", 2*$j+1, " & ",2*$j+2," unequal in length") unless
	  length($seq[2*$j]) == length($seq[2*$j+1]);
# classify and count sites for each sequence pair...
    @seq1 = unpack "(a3)*", $seq[2*$j]; 
    @seq2 = unpack "(a3)*", $seq[2*$j+1];
    %C1 = %C2 = %M = (); # count hashes
    $cdnct = 0;
  CODON:
    foreach $cdn1 (@seq1) {
	$cdn2 = shift @seq2;
	$cdnct++;
	if (($cdn1.$cdn2) =~ /-/) { # strip gaps
	    warn ("dnds: sequence ",2*$j+1," gaps out of frame at codon $cdnct") unless $cdn1 =~ /(---)|(\w\w\w)/;
	    warn ("dnds: sequence ",2*$j+2," gaps out of frame at codon $cdnct") unless $cdn2 =~ /(---)|(\w\w\w)/;
	    next CODON;
	}
	foreach (unpack "(a)*",$sitab->{$cdn1}) {
	    $C1{$_}++;
	}
	foreach (unpack "(a)*",$sitab->{$cdn2}) {
	    $C2{$_}++;
	}
	# deal with %M via $pthtab->{$cdn1.$cdn2}
	#
	if ($cdn1 ne $cdn2) {
	    $pth = $pthtab->{$cdn1.$cdn2};
	    $M{syn} += $pth->{syn};
	    $M{non} += $pth->{non};
	    foreach (qw/0 4 S s x V v C c X w y z/) {
		$M{S}->{$_} += $pth->{S}{$_};
		$M{V}->{$_} += $pth->{V}{$_};
	    }
	}

    }
    # go and calculate dn and ds in various ways here
    # each routine returns statistics in a hash ref
    # stats are different depending on routine
  METHOD:
    for ($opt_M) {
	/N/ && do {$stats = nei($opt_R,\%C1,\%C2,\%M); last METHOD;};
	/L/ && do {$stats = li(\%C1,\%C2,\%M); last METHOD;};
	do { die "dnds: -M option '$opt_M' not recognized" };
    }
    # output stats for current pair here
    1;
    if ($j==0) { #print header
	print join("\t", ("seq1","seq2",keys %{$stats})), "\n";
    }
    $fmt = join("\t", ("%.6g") x scalar keys %{$stats});
    $fmt = "$nm[2*$j]\t$nm[2*$j+1]\t".$fmt;
    printf $fmt, map {$stats->{$_}} keys %{$stats}; 
    print "\n";
}
1;

# nei($tvratio, \%counts1, \%counts2, \%mutncounts)
# calculate dnds stats using Nei-Gojobori method
sub nei {
    my ($R, $C1, $C2, $M) = @_;
    my ($ret,$i,$S,$N,$ts,$ps,$pn);
    
    $ts = $R/($R+1); # prob of transition

    $S = $C1->{4} + $C2->{4} +
	$ts*( $C1->{S} + $C2->{S} ) +
	((2*$ts)/(1+$ts))*( $C1->{s} + $C2->{s} )+
	($C1->{x} + $C2->{x}) +
	(1-$ts)*( $C1->{V} + $C2->{V} ) +
	($C1->{v} + $C2->{v})+
	(1+$ts)*0.5*( $C1->{C} + $C2->{C} ) +
	($C1->{c} + $C2->{c})+	
	(1-$ts)*0.5*( $C1->{X} + $C2->{X} ) +
	0.5*( $C1->{w} + $C2->{w} ) +
	((1-$ts)/(1+$ts))*($C1->{y} + $C2->{y} ) +
	( $C1->{z} + $C2->{z} );
    $S /= 2; # average for both seqs
    $N = $C1->{0} + $C2->{0} +
	$ts*( $C1->{V} + $C2->{V} ) +
	((2*$ts)/(1+$ts))*( $C1->{y} + $C2->{y} )+
	(1-$ts)*( $C1->{S} + $C2->{S} ) +
	(1+$ts)*0.5*( $C1->{X} + $C2->{X} ) +
	(1-$ts)*0.5*( $C1->{C} + $C2->{C} ) +
	0.5*( $C1->{w} + $C2->{w} ) +
	((1-$ts)/(1+$ts))*($C1->{s} + $C2->{s} );
    $N /= 2;
    $ps = $S ? $M->{syn}/$S : 0;
    $pn = $N ? $M->{non}/$N : 0;
    # if $ps, $pn > .75, site is saturated-return -1 in distance
    $ret = {
	'ds' => (($ps < 0.75) ? -(0.75)*log(1- (4/3)*$ps) : -1),
	'dsv' => ((($ps < 0.75) && ($ps >0)) ? ($ps*(1-$ps))/( $S*(1-(4/3)*$ps)*(1-(4/3)*$ps)) : 0),
	'dn' => (($pn < 0.75) ? -(0.75)*log(1- (4/3)*$pn) : -1),
	'dnv' => ((($pn < 0.75) && ($pn >0)) ? ($pn*(1-$pn))/( $N*(1-(4/3)*$pn)*(1-(4/3)*$pn)) : 0)
	};
    return $ret;
}
	
# li(\%counts1, \%counts2, \%mutncounts)
# calculate dnds stats using Li-Wu-Luo method
sub li {
    my ($C1, $C2, $M) = @_;
    my ($ret,$i,@L,@A,@B,@P,@Q,@a,@b,@k,@VA,@VB,@VAB,$ps,$pn);

    # calc 0,2,4fold sites
    $L[0] = 0;
    map {$L[0]+=$C1->{$_}} qw/0/;
    map {$L[0]+=$C2->{$_}} qw/0/;
    $L[0] /= 2;
    $L[2] = 0;
    map {$L[2]+=$C1->{$_}} qw/S s V C X w y x v c z/;
    map {$L[2]+=$C2->{$_}} qw/S s V C X w y x v c z/;
    $L[2] /= 2;
    $L[4] = 0;
    map {$L[4]+=$C1->{$_}} qw/4/;
    map {$L[4]+=$C2->{$_}} qw/4/;
    $L[4] /= 2;

    # calc ts, tv proportions within each class
    $P[0] = 0;
    map {$P[0]+=$M->{S}{$_}} qw/0/;
    $P[0] /= $L[0];

    $P[2] = 0;
    map {$P[2]+=$M->{S}{$_}} qw/S s V C X w y/;
    $P[2] /= $L[2];

    $P[4] = 0;
    map {$P[4]+=$M->{S}{$_}}qw/4 x v c z/;
    $P[4] /= $L[4];

    $Q[0] = 0;
    map {$Q[0]+=$M->{V}{$_}} qw/0/;
    $Q[0] /= $L[0];

    $Q[2] = 0;
    map {$Q[2]+=$M->{V}{$_}} qw/S s V C X w y/;
    $Q[2] /= $L[2];

    $Q[4] = 0;
    map {$Q[4]+=$M->{V}{$_}}qw/4 x v c z/;
    $Q[4] /= $L[4];

    # calc auxiliary values
    foreach $i (0,2,4) {
	$a = $a[$i] = 1/(1-2*$P[$i]-$Q[$i]);
	$b = $b[$i] = 1/(1-2*$Q[$i]);
	$c[$i] = ($a-$b)/2;
	$A[$i] = 0.5*log($a)-0.25*log($b);
	$B[$i] = 0.5*log($b);
	$K[$i] = $A[$i] + $B[$i];
    }
    
    # calc variances
    foreach $i (0,2,4) {
	$VA[$i] = ($a[$i]*$a[$i]*$P[$i] + $c[$i]*$c[$i]*$Q[$i] -
		   ($a[$i]*$P[$i] + $c[$i]*$Q[$i])*($a[$i]*$P[$i] + $c[$i]*$Q[$i]))/$L[$i];
	$VB[$i] = $b[$i]*$b[$i]*$Q[$i]*(1-$Q[$i])/$L[$i];
	$VAB[$i] = ($a[$i]*$a[$i]*$P[$i] + $c[$i]*$c[$i]*$Q[$i] -
		    ($a[$i]*$P[$i] + $k[$i]*$Q[$i])*
		    ($a[$i]*$P[$i] + $k[$i]*$Q[$i]))/$L[$i];
    }
		    

    $ret = {
	'ds' => 3*( $L[2]*$A[2] + $L[4]*$K[4] )/($L[2]+3*$L[4]),
	'dsv' => 9*($L[2]*$L[2]*$VA[2] + $L[4]*$L[4]*$VAB[4])/(($L[2]+3*$L[4])*($L[2]+3*$L[4])),
	'dn' => 3*( $L[2]*$B[2] + $L[0]*$K[0])/(2*$L[2]+3*$L[0]),
	'dnv' => 9*($L[2]*$L[2]*$VB[2] + $L[0]*$L[0]*$VAB[0])/((2*$L[2]+3*$L[0])*(2*$L[2]+3*$L[0])),
	'd4' => $K[4],
	'd4v' => $VAB[4],
	'd0' => $K[0],
	'd0v' => $VAB[0]
    };
    return $ret;
}
