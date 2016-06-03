#!/usr/bin/perl

############## LICENSE ##############
# Copyright (c) 2014-2016 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of cgpCaVEManPostProcessing.
#
# cgpCaVEManPostProcessing is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#    1. The usage of a range of years within a copyright statement contained within
#    this distribution should be interpreted as being equivalent to a list of years
#    including the first and last year specified and all consecutive years between
#    them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
#    2009, 2011-2012’ should be interpreted as being identical to a statement that
#    reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
#    statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
#    identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
#    2009, 2010, 2011, 2012’."
#
############## LICENSE ##############

use strict;
use warnings FATAL=>'all';
use autodie;
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Carp;
use Const::Fast;
use File::Copy qw(move);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

const my $HEAD_FORMAT_1 => "#Unmatched normal panel for CaVEMan";
const my $HEAD_FORMAT_2 => "#PARAMS\tvcfUnmatchedMinMutAlleleCvg=%d\tvcfUnmatchedMinSamplePct=%.3f\tpanelMD5=%s";
const my $BED_LINE_FORMAT => "%s\t%d\t%d\t%.3f\t%s\n";
const my @BASES => ('A','C','G','T');

{
  #Setup options
  my $opts = setup_options();

  my ($md5, $max_col) = verify_data($opts);

  my ($started, $chr_file, $bed_start, $start_file, $score_file, $allele_file) = existing_output($opts, $md5);

  my $conv_prob = $opts->{'p'} / 100; # divide this once, rather than multiply on every record

  open(my $OUT, '>>', $opts->{'out'}) or croak("Error opening output file '".$opts->{'out'}."': $!");
  if($started==0){
    #Print header with params used for generating vcf unmatched normal panel
    print $OUT $HEAD_FORMAT_1,"\n";
    print $OUT sprintf($HEAD_FORMAT_2,$opts->{'c'},$opts->{'p'},$md5 ),"\n" or croak ("Error writing header line to output file.");
  }

  my @base_list;
  my $max_allele = 0;
  my $old_gt = q{};

  my $command = sprintf 'gunzip -c %s', $opts->{'file'};
  my ($pid, $process);
  $pid = open $process, q{-|}, $command or croak 'Could not fork: '.$!;

  my $header_lines = 0;
  while(my $line = <$process>) {
    next if ($line =~ m/^#/);
    chomp ($line);
    my @x = split /\t/, $line;
    my $chr = $x[0];
    my $refb = $x[3];
    my $pos = $x[1];
    if ($started){
      #Skip through VCF until we hit the last line.
      next unless ($chr_file eq $chr && $start_file == $pos);
    }
    my $sample_count = 0;
    my %fail_samples = ('A' => 0,
                        'C' => 0,
                        'G' => 0,
                        'T' => 0,);

    if($#base_list == -1 || $old_gt ne $x[8]){ #Check for differing gt line or not yet set
      my @gt_heads = split /:/, $x[8];
      @base_list = ();
      for (my $j=1;$j<scalar(@gt_heads);$j++){ #Start at one so we skip the genotype
        my $base = $gt_heads[$j];
        $base =~ s/Z//;
        push @base_list, $base;
      }
      $max_allele = scalar @base_list;
      $old_gt = $x[8];
    }

    for my $i(9..$max_col) {
      $sample_count++;
      next if($x[$i] =~m/^[.\-]$/);
      my @gts = split /:/, $x[$i];
      for my $j(1..$max_allele) {
        next if($base_list[$j-1] eq $refb);
        $fail_samples{$base_list[$j-1]}++ if($gts[$j] >= $opts->{'c'});
      } #End of iteration through each base in this genotype
    } #End of iteration through each sample
    my $fail = q{};
    my $score = 0;
    foreach my $base(@BASES){
      next if($base eq $refb);
      my $this_score = ($fail_samples{$base}/$sample_count);
      if($this_score >= $conv_prob){
        $fail .= $base;
        $score = $this_score if($this_score > $score);
      }
    }
    #If this is the last line of an existing file check they match before continuing
    if($started && $start_file == $pos && $chr_file eq $chr){
      croak("Last line in existing file does not match old:'$chr_file,$start_file,$score_file,$allele_file' new:'$chr,$pos,".sprintf("%.3f",$score).",$fail'.") if(sprintf("%.3f",$score) ne $score_file || $fail ne $allele_file);
      print "Found previous line in VCF, starting to analyse again\n";
      next;
    }
    if($. % 1000000 == 0){
      warn "Processed $. lines of VCF.\n";
    }
    printf $OUT $BED_LINE_FORMAT,$chr,$pos-1,$pos,$score,$fail if($fail ne "");

  }#End of parsing each line of data in vcf file
  print "Processed a total of $. lines of VCF.\n";

  close $process or croak "Error closing pipe for $command";

  #End of iteration through each VCF file
  close($OUT);
}

sub existing_output {
  my ($opts, $md5) = @_;
  my $started = 0;
  my $chr_file = q{};
  my $bed_start = 0;
  my $start_file = 0;
  my $score_file = 0;
  my $allele_file = q{};
  my $bad_last_line = 0;
  #If the output file exists,we can try and start
  if(-e $opts->{'out'}){
    $started = 1;
    my $new_file = $opts->{'out'}.'.tmp';
    open (my $NEW, '>', $new_file) or croak ("Error trying to read the ");

    open(my $RD, '<', $opts->{'out'}) or croak ("Error trying to read the ");
    while(<$RD>){
      my $line = $_;
      chomp($line);
      if($line =~ m/^#PARAMS\s*/){
        check_header($opts, $md5, $line);
      }

      next if($line =~ m/^#/);
      my @split = split /\t/, $line;
      if(scalar(@split)!=5){
        #Don't assign partial or bad lines as we can't guarantee the last line was flushed to file fully.
        #If this is the case we mark as to remove the final line.
        $bad_last_line = 1;
        next;
      }else{
        print $NEW $line,"\n" or croak ("Error writing line");
      }
      $chr_file = $split[0];
      $bed_start = $split[1];
      $start_file = $split[2];
      $score_file = $split[3];
      $allele_file = $split[4];
    }
    close($RD);

    close($NEW);
    print "Last line in bed file found: $chr_file, $start_file,$score_file,$allele_file\n";

    if($bad_last_line){ #If we have bad lines we move files from new to the original.
      move($new_file, $opts->{'out'}) or croak("Error trying to move '$new_file' -> '".$opts->{'out'}."': $!");
    }else{
      #Not a bad file so we remove the 'new' file.
      unlink($new_file);
    }

  }#End of checking existing results.
  return ($started, $chr_file, $bed_start, $start_file, $score_file, $allele_file);
}

sub verify_data {
  my ($opts) = @_;
  my $md5;
  my $max_col = 8; # cols is info + genotype cols, but 0 index
  #Get the md5 of samples
  my $Z = new IO::Uncompress::Gunzip $opts->{'file'}, {MultiStream => 1} or die "gunzip failed: $GunzipError\n";
  while(<$Z>){
    my $line = $_;
    next if ($line =~ m/^\s*##/);
    chomp ($line);
    if($line =~ m/^\s*#CHROM/){ #This must be the CHROM line
      my @split_chrom = split /\t/, $line;
      my @samples = ();
      for(my $i=9; $i<scalar(@split_chrom); $i++){
        push @samples, $split_chrom[$i];
      }
      #Sort sample names lexically
      my @sorted = sort @samples;
      $max_col += scalar @sorted;
      my $sample_line = join("\t",@sorted);
      $md5 = md5_hex($sample_line);
      last;
    }
  }
  close($Z);
  return ($md5, $max_col);
}

sub check_header {
  my ($opts, $md5, $line) = @_;
  my @split = split /\s+/, $line;
  foreach my $entry(@split){
    if($entry =~ m/^vcfUnmatchedMinMutAlleleCvg=/){
      my $val = $entry ;
      $val =~ s/vcfUnmatchedMinMutAlleleCvg=//;
      croak("Existing header vcfUnmatchedMinMutAlleleCvg and commandline don't match.") unless($val == $opts->{'c'});
    }
    if($entry =~ m/^vcfUnmatchedMinSamplePct=/){
      my $val = $entry ;
      $val =~ s/vcfUnmatchedMinSamplePct=//;
      croak("Existing header vcfUnmatchedMinSamplePct and commandline don't match.") unless ($val eq (sprintf("%.3f",$opts->{'p'})));
    }
    if($entry =~ m/^panelMD5=/){
      my $val = $entry;
      $val =~ s/panelMD5=//;
      croak("Existing header panelMD5 and commandline don't match.") unless($val eq $md5);
    }
  }
}


sub setup_options {
  my %opts;
  GetOptions(
  				'h|help' => \$opts{'h'},
					'm|man' => \$opts{'m'},
					'f|infile=s' => \$opts{'file'},
					'o|output=s' => \$opts{'out'},
					'c|vcfUnmatchedMinMutAlleleCvg=i' => \$opts{'c'},
					'p|vcfUnmatchedMinSamplePct=f' => \$opts{'p'},
					) or pod2usage(2);

	pod2usage(-verbose => 2) if(defined $opts{'h'});
  pod2usage(-verbose => 1) if(defined $opts{'m'});

  pod2usage(-message => "Option 'f|infile' List of vcf files to convert required.", -verbose => 0) if(!defined $opts{'file'});
  pod2usage(-message => "Option 'o|output' Output location required.", -verbose => 0) if(!defined $opts{'out'});
  pod2usage(-message => "Option 'c|vcfUnmatchedMinMutAlleleCvg' or 'p|vcfUnmatchedMinSamplePct' params required.", -verbose => 0) if(!defined $opts{'c'} || !defined($opts{'p'}));

  return \%opts;
}

__END__

=head1 NAME

convertVCFUnmatchedToBed.pl - Convert unmatched vcf panel to bed file, sort, gzip and tabix.

=head1 SYNOPSIS

convertVCFUnmatchedToBed.pl [options]

  Required parameters:
    -infile                       -f   List of gzipped unmatched vcf files to convert to a single bed file
    -output                       -o   Output bed file
    -vcfUnmatchedMinMutAlleleCvg  -c   Min mutant allele coverage used in flagging;
    -vcfUnmatchedMinSamplePct     -p   Min sample percentage to fail um normal check used in flagging

  Other:
    -help                  -h   Brief help message.
    -man                   -m   Full documentation.

=head1 OPTIONS

=over 8

=item B<-infile>

List of gzipped unmatched vcf files to convert to a single bed file

=item B<-output>

Output file to write bed format to

=item B<-vcfUnmatchedMinMutAlleleCvg>

Parameter used in CaVEMan flagging that this file was built for

=item B<-vcfUnmatchedMinSamplePct>

Parameter used in CaVEMan flagging that this file was built for


=item B<-help>

Prints the help for this script

=item B<-man>

Prints the man page for this script

=back

=head1 DESCRIPTION

B<convertVCFUnmatchedToBed.pl>

=cut



