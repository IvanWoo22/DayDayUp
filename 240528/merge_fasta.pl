#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub READ_FASTA {
    my ($FILE) = @_;
    my %SEQ;
    my $HEADER = '';
    open( my $FH, '<', $FILE ) or die "Could not open file '$FILE' $!";
    while ( my $LINE = <$FH> ) {
        chomp $LINE;
        if ( $LINE =~ /^>(.*)/ ) {
            $HEADER = $1;
        }
        else {
            $SEQ{$HEADER} .= $LINE;
        }
    }
    close $FH;
    return \%SEQ;
}

sub FIND_SEQ {
    my ($SEQ) = @_;
    my @NAMES = keys %$SEQ;
    my %GRAPH;
    for ( my $i = 0 ; $i < @NAMES ; $i++ ) {
        for ( my $j = $i + 1 ; $j < @NAMES ; $j++ ) {
            if (   $SEQ->{ $NAMES[$i] } eq $SEQ->{ $NAMES[$j] }
                || index( $SEQ->{ $NAMES[$i] }, $SEQ->{ $NAMES[$j] } ) == 0
                || index( $SEQ->{ $NAMES[$j] }, $SEQ->{ $NAMES[$i] } ) == 0 )
            {
                push @{ $GRAPH{ $NAMES[$i] } }, $NAMES[$j];
                push @{ $GRAPH{ $NAMES[$j] } }, $NAMES[$i];
                splice( @NAMES, $j, 1 );
                $j--;
            }
        }
    }
    return \%GRAPH;
}

sub DFS {
    my ( $NODE, $VISITED, $GRAPH, $COMPONENT ) = @_;
    $VISITED->{$NODE} = 1;
    push @$COMPONENT, $NODE;
    foreach my $NEIGHBOR ( @{ $GRAPH->{$NODE} } ) {
        if ( !$VISITED->{$NEIGHBOR} ) {
            DFS( $NEIGHBOR, $VISITED, $GRAPH, $COMPONENT );
        }
    }
}

sub FIND_ALL {
    my ($GRAPH) = @_;
    my %VISITED;
    my @COMPONENTS;
    foreach my $NODE ( keys %$GRAPH ) {
        if ( !$VISITED{$NODE} ) {
            my @COMPONENT;
            DFS( $NODE, \%VISITED, $GRAPH, \@COMPONENT );
            push @COMPONENTS, \@COMPONENT;
        }
    }
    return @COMPONENTS;
}

my $fasta_file        = $ARGV[0];
my $sequences         = READ_FASTA($fasta_file);
my $graph             = FIND_SEQ($sequences);
my @components        = FIND_ALL($graph);
my %all_in_components = map { $_ => 1 } map { @$_ } @components;
foreach my $seq_name ( keys %$sequences ) {
    if ( !exists $all_in_components{$seq_name} ) {
        push @components, [$seq_name];
    }
}

my @components_info;
foreach my $component (@components) {
    my $intersection_sequence = $sequences->{ $component->[0] };
    foreach my $seq (@$component) {
        if ( length( $sequences->{$seq} ) < length($intersection_sequence) ) {
            $intersection_sequence = $sequences->{$seq};
        }
    }
    push @components_info,
      {
        component             => $component,
        intersection_sequence => $intersection_sequence,
        length                => length($intersection_sequence)
      };
}

@components_info = sort { $b->{length} <=> $a->{length} } @components_info;
if (@components_info) {
    my $set_number = 1;
    foreach my $info (@components_info) {
        my $set_name              = "set" . $set_number++;
        my $intersection_sequence = $info->{intersection_sequence};
        my $sequence_names        = join( "|", @{ $info->{component} } );
        my $sequences_list =
          join( "|", map { $sequences->{$_} } @{ $info->{component} } );
        my $number_of_sequences = scalar( @{ $info->{component} } );
        print
"$set_name\t$intersection_sequence\t$sequence_names\t$sequences_list\t$number_of_sequences\n";
    }
}
else {
    print "No sequences with identical or contained relationships found.\n";
}

__END__
