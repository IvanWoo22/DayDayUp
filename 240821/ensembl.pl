#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Data::Dumper;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $tree_adaptor = $registry->get_adaptor( "Multi", "compara", "GeneTree" );
my $gene_member_adaptor =
  $registry->get_adaptor( "Multi", "compara", "GeneMember" );
my $genome_db_adaptor =
  $registry->get_adaptor( 'Multi', 'compara', 'GenomeDB' );

my $genome_db = $genome_db_adaptor->fetch_by_name_assembly('homo_sapiens');
my $genes     = $gene_member_adaptor->fetch_all_by_GenomeDB($genome_db);
@$genes = @$genes[ 15400 .. $#$genes ];
foreach my $gene (@$genes) {
    next if ( $gene->biotype_group() ne 'coding' );
    my $tree = $tree_adaptor->fetch_default_for_Member($gene);
    next unless $tree;
    print( $tree->nhx_format("gene_id") );
    print "\n";
}

__END__