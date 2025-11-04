#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use utf8;
use Getopt::Long;

my $start = time;

my ( $hmmquery, $pfam_database, $blastp_database );
my $output_dir = ".";
my @input_files;

GetOptions(
    "hmmquery|h=s"        => \$hmmquery,
    "pfam_database|p=s"   => \$pfam_database,
    "blastp_database|b=s" => \$blastp_database,
    "output_dir|o=s"      => \$output_dir,
    "input_file|i=s@"     => \@input_files,

  )
  or die
"Usage: $0 --hmmquery <hmm_file> --pfam_database NB-ARC.hmm --blastp_database angiosperm_190aa_land_plant_ref_R_genes_classification --output_dir <dir> --input_file <file1> [--input_file <file2> ...]\n";

die "Error: --hmmquery is required\n"                unless defined $hmmquery;
die "Error: At least one --input_file is required\n" unless @input_files;

# Ensure output_dir ends with /
$output_dir =~ s|/$||;    # Remove trailing slash if present
$output_dir .= "/";

my $domainID = $hmmquery;
$domainID =~ s/\..*//;

my $genome_number = 1;

foreach my $object_filename (@input_files) {
    my %hash_object_file = ();
    my $catch_seqio_obj =
      Bio::SeqIO->new( -file => "$object_filename", -format => 'fasta' );
    while ( my $seq_obj = $catch_seqio_obj->next_seq ) {
        my $display_name = $seq_obj->display_name;
        my $seq          = $seq_obj->seq;
        $hash_object_file{"$display_name"} = "$seq";
    }

    my $hmmsearch_output = $output_dir . "/hmmresult.txt";
    my $local_hmmsearch  = system(
"hmmsearch -E 0.0001 --domtblout $hmmsearch_output --domE 0.0001 $hmmquery $object_filename >$output_dir"
          . "/hmmtemp.txt" );

    unless ( open( HMMOUT, "$hmmsearch_output" ) ) {
        print "Cannot open file \"$hmmsearch_output\"\n\n";
        exit;
    }
    my @hmmsearch_hits = <HMMOUT>;
    close HMMOUT;

    my @hmmsearch_hits_name = ();
    foreach my $hmm_line (@hmmsearch_hits) {
        if ( $hmm_line =~ /^#/ ) {
            next;
        }
        else {
            my @split_hmm          = split /\s+/, $hmm_line;
            my $hmmsearch_hit_name = $split_hmm[0];
            push( @hmmsearch_hits_name, $hmmsearch_hit_name );
        }
    }

    my @hmmsearch_hits_seq = ();
    my $count_hmm          = "0";
    foreach my $name_hmmsearch (@hmmsearch_hits_name) {
        if ( exists $hash_object_file{$name_hmmsearch} ) {
            push( @hmmsearch_hits_seq, ">$name_hmmsearch\n" );
            push( @hmmsearch_hits_seq, "$hash_object_file{$name_hmmsearch}\n" );
        }
        else {
            next;
        }
        $count_hmm += 1;
    }
    print "hmmsearch found $count_hmm hits in $object_filename\n";

    if ( $count_hmm == 0 ) {
        next;
    }
    else {
        my $out_hmm_seq = $output_dir . "/hmmsearch_out.fas";
        open OUT1, ">$out_hmm_seq";
        print OUT1 @hmmsearch_hits_seq;
        close OUT1;

        my $setupDB = system("makeblastdb -in $object_filename -dbtype prot");
        my $blast_result_output = $output_dir . "/blastresult.txt";
        my $local_blast         = system(
"blastp -db $object_filename -query $out_hmm_seq -out $blast_result_output -evalue 0.0001 -outfmt 6"
        );

        unless ( open( BLASTRESULTS, "$blast_result_output" ) ) {
            print "Cannot open file \"$blast_result_output\"\n\n";
            exit;
        }
        my @Blast_hits = <BLASTRESULTS>;
        close BLASTRESULTS;

        my @Blast_hits_name = ();
        foreach (@Blast_hits) {
            my @split_blast    = split /\t/, $_;
            my $hit_name_blast = $split_blast[1];
            push( @Blast_hits_name, $hit_name_blast );
        }

        my @Blast_hits_seq = ();
        my $count_blast    = "0";
        my %unique         = ();
        foreach my $item (@Blast_hits_name) {
            next unless defined $item;    # 跳过未定义的元素
            $unique{$item}++;
        }
        my @unique_blast_hits_name = keys %unique;

        foreach my $name_blast (@unique_blast_hits_name) {
            if ( exists $hash_object_file{$name_blast} ) {
                push( @Blast_hits_seq, ">$name_blast\n" );
                push( @Blast_hits_seq, "$hash_object_file{$name_blast}\n" );
            }
            else {
                next;
            }
            $count_blast += 1;
        }
        print "Blast found $count_blast hits in $object_filename\n";

        my $out_blast_seq = $output_dir . "/blast_seq_out.fas";
        open OUT2, ">$out_blast_seq";
        print OUT2 @Blast_hits_seq;
        close OUT2;

        my $hmmscan_output = $output_dir . "/hmmscan_output.txt";
        my $local_hmmscan  = system(
"hmmscan -E 1e-4 --incE 1e-4 --domtblout $hmmscan_output  $pfam_database $out_blast_seq >$output_dir"
              . "/hmmtemp.txt" );
        unless ( open( HMMSCANOUT, "$hmmscan_output" ) ) {
            print "Cannot open file \"$hmmscan_output\"\n\n";
            exit;
        }
        my @Pfam_hits = <HMMSCANOUT>;
        close HMMSCANOUT;

        my @hmmscan_hits_name = ();
        foreach my $pfam_hit_name (@Pfam_hits) {
            my @split_hmmscan = split /\s+/, $pfam_hit_name;
            my $domain        = $split_hmmscan[0];
            my $hit_name      = $split_hmmscan[3];
            my $evalue        = $split_hmmscan[12];
            if ( $domain =~ /\Q$domainID\E/ && $evalue <= 0.0001 ) {
                push( @hmmscan_hits_name, $hit_name );
            }
        }

        my %unique_name = ();
        foreach my $single (@hmmscan_hits_name) {
            $unique_name{$single}++;
        }
        my @unique_hmmscan_hits_name = keys %unique_name;

        my %hash_blast_out = ();
        my $catch_seqio_obj_new =
          Bio::SeqIO->new( -file => "$object_filename", -format => 'fasta' );
        while ( my $seq_obj = $catch_seqio_obj_new->next_seq ) {
            my $display_name = $seq_obj->display_name;
            my $seq          = $seq_obj->seq;
            $hash_blast_out{"$display_name"} = "$seq";
        }

        my @hmmscan_hits_seq = ();
        my $count_hmmscan    = "0";
        foreach my $uni_name_hmmscan (@unique_hmmscan_hits_name) {
            if ( exists $hash_blast_out{$uni_name_hmmscan} ) {
                push( @hmmscan_hits_seq, ">$uni_name_hmmscan\n" );
                push( @hmmscan_hits_seq,
                    "$hash_blast_out{$uni_name_hmmscan}\n" );
            }
            else {
                next;
            }
            $count_hmmscan += 1;
        }
        print "Local Pfam found $count_hmmscan hits in $object_filename\n";

        my $out_hmmscan_seq = $output_dir . "/hmmscan_seq_out.fas";
        open OUT3, ">$out_hmmscan_seq";
        print OUT3 @hmmscan_hits_seq;
        close OUT3;

        #classification
        my $blast_m8_file   = $output_dir . "/1st_blast_classification.txt";
        my $local_blast_1st = system(
"blastp -db $blastp_database -query $out_hmmscan_seq -out $blast_m8_file -evalue 0.001 -max_target_seqs 1 -outfmt 6"
        );

        #classification_over
        print "The $genome_number genome is completed\n";
        $genome_number += 1;
        print "\t\tThe script is running for: ", time - $start, " seconds\n";
    }
}
