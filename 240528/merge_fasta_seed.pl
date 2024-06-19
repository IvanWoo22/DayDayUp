#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

# Function to read a FASTA file and return a hash reference of sequences
sub READ_FASTA {
    my ($FILE) = @_;
    my %SEQUENCES;
    my $HEADER = '';
    open my $FH, '<', $FILE;
    while ( my $LINE = <$FH> ) {
        chomp $LINE;
        if ( $LINE =~ /^>(.*)/ ) {
            $HEADER = $1;
        }
        else {
            $SEQUENCES{$HEADER} .= $LINE;
        }
    }
    close $FH;
    return \%SEQUENCES;
}

# Function to group sequences by their first 15 characters
sub GROUP_BY_PREFIX {
    my ($SEQUENCES) = @_;
    my %GROUPS;
    for my $NAME ( keys %$SEQUENCES ) {
        my $PREFIX = substr( $SEQUENCES->{$NAME}, 0, 15 );
        push @{ $GROUPS{$PREFIX} }, $NAME;
    }
    return \%GROUPS;
}

# Function to find sequences that are identical or one is a prefix of the other
sub FIND_SEQUENCES {
    my ($SEQUENCES) = @_;
    my @NAMES = keys %$SEQUENCES;
    my %GRAPH;
    for ( my $I = 0 ; $I < @NAMES ; $I++ ) {
        for ( my $J = $I + 1 ; $J < @NAMES ; $J++ ) {
            if ( $SEQUENCES->{ $NAMES[$I] } eq $SEQUENCES->{ $NAMES[$J] }
                || index( $SEQUENCES->{ $NAMES[$I] },
                    $SEQUENCES->{ $NAMES[$J] } ) == 0
                || index( $SEQUENCES->{ $NAMES[$J] },
                    $SEQUENCES->{ $NAMES[$I] } ) == 0 )
            {
                push @{ $GRAPH{ $NAMES[$I] } }, $NAMES[$J];
                push @{ $GRAPH{ $NAMES[$J] } }, $NAMES[$I];
                splice @NAMES, $J, 1;
                $J--;
            }
        }
    }
    return \%GRAPH;
}

# Depth-First Search to explore graph components
sub DFS {
    my ( $NODE, $VISITED, $GRAPH, $COMPONENT ) = @_;
    $VISITED->{$NODE} = 1;
    push @$COMPONENT, $NODE;
    for my $NEIGHBOR ( @{ $GRAPH->{$NODE} } ) {
        DFS( $NEIGHBOR, $VISITED, $GRAPH, $COMPONENT )
          unless $VISITED->{$NEIGHBOR};
    }
}

# Find all connected components in the graph
sub FIND_ALL_COMPONENTS {
    my ($GRAPH) = @_;
    my %VISITED;
    my @COMPONENTS;
    for my $NODE ( keys %$GRAPH ) {
        unless ( $VISITED{$NODE} ) {
            my @COMPONENT;
            DFS( $NODE, \%VISITED, $GRAPH, \@COMPONENT );
            push @COMPONENTS, \@COMPONENT;
        }
    }
    return @COMPONENTS;
}

# Main script logic
my $fasta_file = $ARGV[0];
my $sequences  = READ_FASTA($fasta_file);
my $groups     = GROUP_BY_PREFIX($sequences);

my @all_components;
for my $prefix ( keys %$groups ) {
    my %group_sequences =
      map { $_ => $sequences->{$_} } @{ $groups->{$prefix} };
    my $group_graph = FIND_SEQUENCES( \%group_sequences );
    my @components  = FIND_ALL_COMPONENTS($group_graph);
    push @all_components, @components;
}

my %all_in_components = map { $_ => 1 } map { @$_ } @all_components;
for my $seq_name ( keys %$sequences ) {
    push @all_components, [$seq_name] unless $all_in_components{$seq_name};
}

# Collect information about components
my @components_info;
for my $component (@all_components) {
    my $intersection_sequence = $sequences->{ $component->[0] };
    for my $seq (@$component) {
        $intersection_sequence = $sequences->{$seq}
          if length( $sequences->{$seq} ) < length($intersection_sequence);
    }
    push @components_info,
      {
        COMPONENT             => $component,
        INTERSECTION_SEQUENCE => $intersection_sequence,
        LENGTH                => length($intersection_sequence),
      };
}

# Sort and print components info
@components_info = sort { $b->{LENGTH} <=> $a->{LENGTH} } @components_info;
if (@components_info) {
    my $set_number = 1;
    for my $info (@components_info) {
        my $set_name       = sprintf( "set%08d", $set_number++ );
        my $sequence_names = join( "|", @{ $info->{COMPONENT} } );
        my $sequences_list =
          join( "|", map { $sequences->{$_} } @{ $info->{COMPONENT} } );
        print
"$set_name\t$info->{INTERSECTION_SEQUENCE}\t$sequence_names\t$sequences_list\t",
          scalar( @{ $info->{COMPONENT} } ), "\n";
    }
}
else {
    print "No sequences with identical or contained relationships found.\n";
}

__END__