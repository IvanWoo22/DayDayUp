#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# 默认参数
my $MAX_GAP        = 200000;
my $MAX_NON_TARGET = 8;
my $H2H_THR        = 15000;

my $allfile;                    # BED file (sorted)
my $targetfile;                 # target gene list
my $output = "clusters.tsv";    # 输出文件

my $help = 0;
my $man  = 0;

GetOptions(
    "bed|b=s"     => \$allfile,
    "target|t=s"  => \$targetfile,
    "output|o=s"  => \$output,
    "max_gap|g=i" => \$MAX_GAP,
    "max_non|n=i" => \$MAX_NON_TARGET,
    "h2h_thr|h=i" => \$H2H_THR,
    "help|?"      => \$help,
    "man"         => \$man,
) or pod2usage(2);

pod2usage(1)                              if $help != 0;
pod2usage( -exitval => 0, -verbose => 2 ) if $man != 0;

pod2usage("Error: --bed file is required\n")    unless $allfile ne '';
pod2usage("Error: --target file is required\n") unless $targetfile ne '';

# 确保输出目录存在
if ( $output =~ m{(.*)/} ) {
    my $dir = $1;
    unless ( -d $dir ) {
        mkdir $dir or die "Cannot create output directory '$dir': $!\n";
    }
}

# read target list
open my $TF, '<', $targetfile
  or die "Cannot open target file '$targetfile': $!\n";
my %is_target;
while (<$TF>) { chomp; next if /^\s*$/ || /^#/; $is_target{$_} = 1 }
close $TF;

open my $AF, '<', $allfile or die "Cannot open BED file '$allfile': $!\n";
open my $OUTC, '>', $output
  or die "Cannot write to output file '$output': $!\n";

print $OUTC join( "\t",
    qw(ClusterID Chr Start End Size Target_genes Category H2H_flag H2H_distance)
) . "\n";

my $curr_chr              = "";
my $in_cluster            = 0;
my $cluster_id            = 0;
my @cluster_targets_info  = ();    # store "gene|start|end|strand"
my $last_target_end       = 0;
my $non_target_since_last = 0;
my $cluster_start         = 0;

sub emit_cluster {
    my ( $cid, $chr, $cstart, $cend, $targets_ref ) = @_;
    my @targets = @$targets_ref;
    my $size    = scalar @targets;
    my $category =
      $size == 1 ? "Singleton" : ( $size == 2 ? "Paired" : "Cluster" );

    # compute H2H info only for size==2
    my $h2h_flag = "";
    my $h2h_dist = "";
    if ( $size == 2 ) {
        my ( undef, $s1, $e1, $st1 ) = split( /\|/, $targets[0] );
        my ( undef, $s2, $e2, $st2 ) = split( /\|/, $targets[1] );

# TSS: for '+' use start, for '-' use end (GFF->bed start is 0based; using these coords consistently)
        my $tss1 = ( $st1 eq '+' ) ? $s1 : $e1;
        my $tss2 = ( $st2 eq '+' ) ? $s2 : $e2;
        $h2h_dist = abs( $tss1 - $tss2 );
        if ( $st1 ne $st2 && $h2h_dist <= $H2H_THR ) {
            $h2h_flag = "H2H_tight";
        }
        elsif ( $st1 ne $st2 ) {
            $h2h_flag = "H2H_loose";
        }
        else {
            $h2h_flag = "not_H2H";
        }
    }
    print $OUTC join( "\t",
        $cid, $chr, $cstart, $cend, $size,
        join( ",", map { ( split( /\|/, $_ ) )[0] } @targets ),
        $category, $h2h_flag, $h2h_dist )
      . "\n";
}

while (<$AF>) {
    chomp;
    next if /^\s*$/ || /^#/;
    my ( $chr, $start, $end, $gene, undef, $strand ) = split(/\t/);

    # ensure strand present
    $strand = ( defined $strand && $strand ne "" ) ? $strand : ".";
    if ( $chr ne $curr_chr ) {

        # emit pending cluster if any
        if ( $in_cluster != 0 ) {
            $cluster_id++;
            emit_cluster( $cluster_id, $curr_chr, $cluster_start,
                $last_target_end, \@cluster_targets_info );
        }

        # reset
        $curr_chr              = $chr;
        $in_cluster            = 0;
        @cluster_targets_info  = ();
        $non_target_since_last = 0;
        $last_target_end       = 0;
    }
    my $is_t = $is_target{$gene} ? 1 : 0;
    if ( $is_t != 0 ) {
        if ( $in_cluster == 0 ) {
            $in_cluster           = 1;
            $cluster_start        = $start;
            $last_target_end      = $end;
            @cluster_targets_info = ();
            push @cluster_targets_info,
              join( "|", $gene, $start, $end, $strand );
            $non_target_since_last = 0;
        }
        else {
            my $dist = $start - $last_target_end;
            if (   $dist <= $MAX_GAP
                && $non_target_since_last <= $MAX_NON_TARGET )
            {
                push @cluster_targets_info,
                  join( "|", $gene, $start, $end, $strand );
                $last_target_end       = $end;
                $non_target_since_last = 0;
            }
            else {
                # close previous
                $cluster_id++;
                emit_cluster( $cluster_id, $curr_chr, $cluster_start,
                    $last_target_end, \@cluster_targets_info );

                # start new
                @cluster_targets_info = ();
                push @cluster_targets_info,
                  join( "|", $gene, $start, $end, $strand );
                $cluster_start         = $start;
                $last_target_end       = $end;
                $non_target_since_last = 0;
                $in_cluster            = 1;
            }
        }
    }
    else {
        if ( $in_cluster != 0 ) {
            $non_target_since_last++;
        }
    }
}

# EOF: emit last cluster
if ( $in_cluster && @cluster_targets_info ) {
    $cluster_id++;
    emit_cluster( $cluster_id, $curr_chr, $cluster_start, $last_target_end,
        \@cluster_targets_info );
}

close $AF;
close $OUTC;

print "Cluster file generated: '$output'\n";
print "Parameters used:\n";
print "  MAX_GAP         = $MAX_GAP\n";
print "  MAX_NON_TARGET  = $MAX_NON_TARGET\n";
print "  H2H_THRESHOLD   = $H2H_THR\n";

__END__

=head1 NAME

make_clusters.pl - Identify gene clusters from sorted BED and target list

=head1 SYNOPSIS

  perl make_clusters.pl --bed all_genes.sorted.bed --target target_list.txt [options]

=head1 OPTIONS

=over 8

=item B<--bed, -b> <file>

Sorted BED file (6-column: chr, start, end, gene, ., strand). Required.

=item B<--target, -t> <file>

List of target gene names (one per line). Required.

=item B<--output, -o> <file>

Output TSV file path. Default: clusters.tsv

=item B<--max_gap, -g> <int>

Maximum gap (bp) between target genes to merge. Default: 200000

=item B<--max_non, -n> <int>

Maximum non-target genes allowed between targets. Default: 8

=item B<--h2h_thr, -h> <int>

Head-to-head distance threshold for "tight" H2H. Default: 15000

=item B<--help, -?>

Show this help message.

=item B<--man>

Show full manual.

=back

=head1 OUTPUT COLUMNS

  ClusterID, Chr, Start, End, Size, Target_genes, Category, H2H_flag, H2H_distance

=head1 AUTHOR

Wrote by Ivan and improved by Grok.

=cut
