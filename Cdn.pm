package Cdn;

use strict;
use warnings;

BEGIN {
    use Exporter ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);        # set the version for version checking
    $VERSION     = 1.00;
    @ISA         = qw(Exporter);
    @EXPORT      = qw(&path_table &count_sites &cdnpath &read_codon_table);
    %EXPORT_TAGS = ( );
    @EXPORT_OK   = ( );
}
our @EXPORT_OK;

#  count diffs bwn codons, generate paths
#  create as lists single nt change codon paths

# path_table( \%cdntab, \%sitab )
sub path_table {
    my ($cdntab, $sitab) = @_;
    my (%pthtab, @cdn, $cdn, $Tcdn, $npaths, $path, @p,$c,$Rx, $Yx, $Ix,$i,$j);
    # counts
    my ($S, $N, %S, %V);
    $Rx = "A"^"G";
    $Yx = "C"^"T";
    $Ix = "A"^"A";
    
    @cdn = reverse sort keys %{$cdntab};

    while ($Tcdn = shift @cdn) {
	next if $cdntab->{$Tcdn} eq "*"; # skip stop codons
      CODON:
	foreach $cdn (@cdn) {
	    next CODON if $cdntab->{$cdn} eq "*"; # skip stop codons
	    ($npaths, $path) = cdnpath($cdn, $Tcdn, $cdntab); 
            # need to count numbers of changes in paths
	    $S = $N = 0;
	    %S = %V = ();
	    for ($i=0; $i<$npaths; $i++) {
		@p = split / /, $path->[$i];
		for ($c=shift @p; @p; $c=shift@p) {
		    if ($cdntab->{$c} eq $cdntab->{$p[0]}) {
			$S++;
		    }
		    else {
			$N++;
		    }
		    $_ = $c^$p[0];
		    /[^$Ix]/g; $j = pos()-1;# find the change
		      SWITCH: for (substr($_,$j,1)) {
			  # is transition
			  ($_ eq $Rx) || ($_ eq $Yx) &&
			      do { $S{substr($sitab->{$c},$j,1)}++; 
				   last SWITCH; };
			  # else is transversion
			  do { $V{substr($sitab->{$c},$j,1)}++; 
			       last SWITCH; };
		      };
		}
	    }
	    # average numbers over paths...
	    $S /= $npaths; $N /= $npaths;
	    map {$S{$_} /= $npaths} keys %S;
	    map {$V{$_} /= $npaths} keys %V;
	    # table entry
	    $pthtab{$cdn.$Tcdn} = 
	    $pthtab{$Tcdn.$cdn} = { 'syn' => $S,
				    'non' => $N,
				    'S' => {%S},
				    'V' => {%V} };
	    
	}
    }
    return {%pthtab}; # ref to copy of pthtab
}

# cdnpath($cdn, $Tcdn, \%cdntab)
sub cdnpath {
    my ($cdn, $Tcdn, $cdntab) = @_;
    my (@path, $npaths, @cpath);
    cdnpath_guts($cdn, $Tcdn, $cdntab, \@path,\$npaths,\@cpath);
    return ($npaths, [@path]);
}


# cdnpath_guts( $cdn, $Tcdn, \%cdntab, \@path, \$npaths, \@cpath)
sub cdnpath_guts {
# $path[$i] contains codons separated by spaces, from initial to final
#  codon in single-nt change steps
# $npaths contains number of paths
# call with initial codon
# 
    my ($cdn, $Tcdn,  $cdntab, $path, $npaths, $cpath) = @_;
    my ($i, $Mcdn);
    if ($cdntab->{$cdn} eq "*") {  # elim paths containing stops
	pop @$cpath;
	return;
    }
    if ($cdn eq $Tcdn) {
	$path->[$$npaths++] = join(" ",@$cpath,$cdn);
	pop @$cpath;
	return;
    }
    for ($i=0; $i < 3; $i++) {
	if (substr($cdn,$i,1) ne substr($Tcdn,$i,1)) {
	    push @$cpath, $cdn;
	    $Mcdn = $cdn;
	    substr($Mcdn,$i,1) = substr($Tcdn,$i,1);
	    cdnpath_guts($Mcdn,$Tcdn,$cdntab,$path,$npaths,$cpath);
	}

    }
    pop @$cpath;
    return;
}


# count_sites( $sequence, \%cdntab )
# count 0-, 2-, and 4-fold sites for a given (sequence of) codons
sub count_sites {
    my ($seq, $cdntab) = @_;
    my (@seq, $cdn, $b, $ts, @tv, $i, $sit,$mcdn,$SS, $SV, $XS, $XV);
    my $ret = "";

    $seq = uc $seq;
    $seq =~ tr/U/T/;
    die "count_sites: input not nucleic" if $seq =~ /[^ATCG-]/;
    @seq = unpack "(a3)*", $seq; # get into codons
    foreach $cdn (@seq) {
	$sit = "---";
	foreach $i (0..2) {
	    # get transitions and transversions for this base
	    $_ = substr($cdn,$i,1);
	    if (/[AG]/) {
		"AG" =~ /([^$_])/;
		$ts = $1;
		@tv = ("C", "T");
	    }
	    else {
		"CT" =~ /([^$_])/;
		$ts = $1;
		@tv = ("A", "G");
	    }
	    # synonymous ts and tv
	    $SS = $SV = $XS = $XV = 0;
	    $mcdn = $cdn;
	    substr($mcdn, $i, 1) = $ts;
	    if ($cdntab->{$mcdn} eq "*") {
		$XS++;
	    }
	    else {
		$SS += ($cdntab->{$cdn} eq $cdntab->{$mcdn});
	    }
	    substr($mcdn, $i, 1) = $tv[0];
	    if ($cdntab->{$mcdn} eq "*") {
		$XV++;
	    }
	    else {
		$SV += ($cdntab->{$cdn} eq $cdntab->{$mcdn});
	    }
	    substr($mcdn, $i, 1) = $tv[1];
	    if ($cdntab->{$mcdn} eq "*") {
		$XV++;
	    }
	    else {
		$SV += ($cdntab->{$cdn} eq $cdntab->{$mcdn});
	    }

	  SWITCH: {
	      ($SS == 0) && do { 
		  ($SV == 0) && do {substr($sit,$i,1) = "0"; last SWITCH; };
		  ($SV == 1) && do {
		      ($XS == 0 && $XV == 0) && do {substr($sit,$i,1) = "X"; last SWITCH; };
		      ($XS == 1 && $XV == 0) && do {substr($sit,$i,1) = "w"; last SWITCH; };
		      ($XS == 0 && $XV == 1) && do {substr($sit,$i,1) = "y"; last SWITCH; };
		      ($XS == 1 && $XV == 1) && do {substr($sit,$i,1) = "z"; last SWITCH; };
		  };
		  ($SV == 2) && do { 
		      ($XS == 0) && do {substr($sit,$i,1) = "V"; last SWITCH; };
		      ($XS == 1) && do {substr($sit,$i,1) = "v"; last SWITCH; };
		  };
	      };
	      ($SS == 1) && do { 
		  ($SV == 0) && do {
		      ($XV == 0) && do {substr($sit,$i,1) = "S"; last SWITCH; };
		      ($XV == 1) && do {substr($sit,$i,1) = "s"; last SWITCH; };
		      ($XV == 2) && do {substr($sit,$i,1) = "x"; last SWITCH; };
		  };
		  ($SV == 1) && do {
		      ($XV == 0) && do {substr($sit,$i,1) = "C"; last SWITCH; };
		      ($XV == 1) && do {substr($sit,$i,1) = "c"; last SWITCH; };
		  };
		  ($SV == 2) && do { substr($sit,$i,1) = "4"; last SWITCH; };
	      };
	  }
	}
    $ret .= $sit;
    }
    return $ret;
}

    
# read_codon_table( $file )
#read a codon table file in EMBOSS distribution format.
#return a reference to the hash
sub read_codon_table {
    my (@file,$file,@a, $i, %h, %tbl);
    $file = shift;
    open TB, $file or die "read_codon_table:$!";
    @file = <TB>;
    close TB;
    @file = grep /=/, @file;
    chomp @file;
    foreach (@file) {
	@a = split (/\s+=\s+/);
	$h{$a[0]} = [split //, $a[1]];
    }
    foreach $i (0..@{$h{AAs}}-1) {
	$tbl{$h{Base1}[$i].$h{Base2}[$i].$h{Base3}[$i]} = $h{AAs}[$i];
    }
    return {%tbl};
}
	
1;
