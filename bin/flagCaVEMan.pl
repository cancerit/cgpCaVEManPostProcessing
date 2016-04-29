#! /usr/bin/perl

##########LICENCE##########
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
##########LICENCE##########

use strict;
use warnings FATAL=>'all';
use autodie;
use Carp qw(croak cluck);
use Getopt::Long qw(:config pass_through);
use Pod::Usage;

use Bio::DB::HTS::Tabix;

use Config::IniFiles;
use File::ShareDir qw(dist_dir);
use Const::Fast qw(const);
use Try::Tiny qw(try catch);

use FindBin qw($Bin);
use lib "$FindBin::Bin/../lib";

use Sanger::CGP::CavemanPostProcessing;
use Sanger::CGP::CavemanPostProcessing::Config;
use Sanger::CGP::CavemanPostProcessing::CaveVCF;
use Sanger::CGP::CavemanPostProcessing::Flagger;

const my @PERMITTED_AMPLICON_TYPES => qw/followup amplicon/;
const my @PERMITTED_RNA_TYPES => qw/rnaseq rna_seq rna-seq rna/;
const my @PERMITTED_WXS_TYPES => qw/pulldown exome exomic wxs targeted targetted/;
const my @PERMITTED_WGS_TYPES => qw/genome wgs genomic/;
const my @PERMITTED_SEQ_TYPES => (@PERMITTED_WGS_TYPES, @PERMITTED_WXS_TYPES, @PERMITTED_AMPLICON_TYPES, @PERMITTED_RNA_TYPES);

const my $WGS_TYPE => 'WGS';
const my $WXS_TYPE => 'WXS';
const my $AMPLICON_TYPE => 'AMPLICON';
const my $RNASEQ_TYPE => 'RNASEQ';
const my $DEFAULT_LINE_CACHE => 2000;
const my $FLAG_TO_VCF_CONFIG => '%s/flag.to.vcf.convert.ini';
const my $FLAG_CONFIG => '%s/human/flag.vcf.config.ini';
const my $FLAG_TYPE => 'Flag';

const my $VCF_UM_FLAG_NAME => 'unmatchedNormalVcfFlag';
const my $UNMATCHED_FORMAT_VCF => 'VCF';
const my $UNMATCHED_FORMAT_BED => 'BED';

const my $DEFAULT_FLAG_MODULE => 'Sanger::CGP::CavemanPostProcessing::Flagger';

my $intersectFlagStore;
my $umformat;
my $opts = option_builder();
validateInput($opts);
main($opts);

sub main {
  my ($opts) = @_;
  my @lineCache = ();
  #Get flagList and info from config files
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  $cfg->init_config($opts);
  my $flagger = $opts->{'mod'}->new();
  $flagger->init($opts->{'m'},$opts->{'n'},$cfg);
  #Open VCF output file
  warn "Opening input and output VCF '".$opts->{'f'}."' & '".$opts->{'o'}."'" if($opts->{'loud'});
  my $caveVCF = Sanger::CGP::CavemanPostProcessing::CaveVCF->new();
  $caveVCF->init_vcfs($opts->{'f'},$opts->{'o'});
  warn "Writing output VCF header \n" if($opts->{'loud'});
  $caveVCF->output_header($opts,$cfg);
  warn "Performing initial intersects with bedfiles\n" if($opts->{'loud'});
  my $has_um_vcf = 0;
  my $umNormVcf;
  my $vcf_flag;
  foreach my $flag(@{$cfg->flags}){
    if($flag->name eq $VCF_UM_FLAG_NAME){
      $has_um_vcf = 1;
      $vcf_flag = $flag;
      next;
    }
    next unless($flag->is_intersect);
    #Run intersects for each of the potential flags, will skip if there's no file.
		warn "Performing ".$flag->name." intersect\n" if($opts->{'loud'});
		getIntersectMatches($opts->{'f'},$flag->intersect_file,$flag->id);
  }
  if($has_um_vcf){
    warn "Unmatched panel being initialised\n" if($opts->{'loud'});
    $umNormVcf = buildUnmatchedVCFFileListFromReference($opts->{'umv'},$opts->{'ref'},$cfg,$vcf_flag);
  }
	warn "Starting flagging\n" if($opts->{'loud'});

  #Iterate through each VCF variant line
	while (my $x=$caveVCF->input_vcf->next_data_array()){
	  $$x[6]='.'; #This resets all flags.
    #Reset info flags
    foreach my $flag(@{$cfg->flags}){
      if($flag->is_info){
        $$x[7]=$caveVCF->input_vcf->add_info_field($$x[7],$flag->id=>undef);
      }
    }
    my ($results,$info) = flagVCFLine($x,$cfg,$flagger,$umNormVcf);
    #Add the relevant filters or PASS to the filter section.
    $$x[6]=$caveVCF->input_vcf->add_filter($$x[6],%$results);
    $$x[7]=$caveVCF->input_vcf->add_info_field($$x[7],%$info);
    #Validate this line and append this line to the output file;
    push(@lineCache,$caveVCF->input_vcf->format_line($x));
    #Empty line cache if it is at/over the limit set
    if(scalar(@lineCache)>=$opts->{'l'}){
      $caveVCF->output_vcf_lines(\@lineCache);
      @lineCache = ();
    }#End of emptying line cache
	}#End of iteration through each variant

	$caveVCF->output_vcf_lines(\@lineCache);

  #Close vcf output file
  $caveVCF->close_vcfs();
  warn "Flagging Done\n" if($opts->{'loud'});

}

sub flagVCFLine{
  my ($x,$cfg,$flagger,$umNormVcf) = @_;
  $flagger->set_position($$x[0],$$x[1],$$x[3],$$x[4]);
  my $results;
  my $info;
  foreach my $flag(@{$cfg->flags}){
    my $name = $flag->name;
    my $res;
    if($flag->is_intersect){
      #Intersect already done
      my $coord = $$x[0].':'.$$x[1];
      if(exists($intersectFlagStore->{$coord}->{$flag->id})){
        $res = 1;
      }
      #croak("intersect not implemented yet");
    }elsif($flag->name eq $VCF_UM_FLAG_NAME){
      my $match = getUnmatchedVCFIntersectMatch($$x[0],$$x[1],$umNormVcf->{$$x[0]});
      $res = 0;
      if(defined($match)){
        #Unmatched VCF normal flag
        if($umformat eq $UNMATCHED_FORMAT_VCF){
          $res = $flagger->getVCFUnmatchedFlag($match);
        }elsif($umformat eq $UNMATCHED_FORMAT_BED){ #Bed format
          $res = $flagger->getBEDUnmatchedFlag($match);
        }else{
          croak("Unknown format '$umformat' for UMVcf");
        }
      }#End of if we have a matching unmatched panel entry
    }else{#A non intersect filter to run
      $res = $flagger->$name();
    }
    #Add the result to the appropriate store
    if($flag->is_info){
      $res='' if(uc($flag->type) eq uc($FLAG_TYPE) && defined $res && $res==1);
      $info->{$flag->id} = $res;
    }else{
      $results->{$flag->id} = $res;
    }
  }
  return ($results,$info);
}

sub getUnmatchedVCFIntersectMatch{
	my ($chr,$pos,$tabix) = @_;
	#Only check for positions we've not already sorted.
	return undef unless defined($tabix);
  # new tabix is 1-based for both coordinates
  my $iter = $tabix->query(sprintf '%s:%d-%d', $chr,$pos,$pos);
  my $line = $iter->next; # undef if none
  return $line;
}

sub check_exists_remote{
  my ($fileloc) = @_;
  return 0 if(!defined($fileloc));
  if($fileloc =~ /^(http|ftp)/){
    if (head($fileloc)) {
      return 1;
    }else{
      return 0;
    }
  }else{
    return 1 if(-e $fileloc);
  }
  return 0;
}

sub buildUnmatchedVCFForSingleBed{
  my ($bedloc,$cfg,$refFai) = @_;
  my $fileList;
  my $REF;
  #Check here that the header matches the params set in the configs.
  croak("Unmatched panel parameters in bed provided don't match those in the config file.") if(_checkConfigAndBedUnmatched($bedloc,$cfg)!=1);
  #Check the tabix file exists.
  my $umBedTabix = $bedloc.".tbi";
  croak("Unmatched bedfile $umBedTabix") if (! -e $umBedTabix);
  open($REF, '<', $refFai) or croak("Error opening reference index $refFai: $!");
    while(<$REF>){
      my $line = $_;
      next if($line =~ m/^\s*#/);
      chomp($line);
      my ($chr,undef) = split(/\t/,$line);
      $fileList->{$chr} = Bio::DB::HTS::Tabix->new(filename => $bedloc);
    }
  close($REF);
  return $fileList;
}

sub buildUnmatchedVCFForMultiVCF{
  my ($bedloc,$umLoc,$cfg,$refFai,$umformat) = @_;
  my $fileList;
  my $count=0;
  my $REF;
  open($REF, '<', $refFai) or croak("Error opening reference index $refFai: $!");
    while(<$REF>){
      my $line = $_;
      next if($line =~ m/^\s*#/);
      chomp($line);
      my ($chr,undef) = split(/\t/,$line);
      $count++;
      #The files will be named according to line index rather than contig name.
      my $umVcfTabix;
      if($umformat eq $UNMATCHED_FORMAT_VCF){
        my $fileName = $umLoc."/unmatchedNormal.".$count.".vcf.gz";
        if (!check_exists_remote($fileName)){
          croak("Couldn't find vcf file '$fileName'.");
        }
        if($fileName =~ /^(http|ftp)/) {
          $umVcfTabix = Bio::DB::HTS::Tabix->new(filename => $fileName);
        }
        else {
          #If file exists, add it to the list, otherwise, fail.
          if(! -e $fileName){
            warn("Couldn't find the unmatched VCF normal file $fileName corresponding to contig $chr in given location.\nAny mutations on this contig will NOT be flagged with the unmatched VCF normal flag.");
            next;
          }
          my $idx = $fileName.".tbi";
          croak ("Tabix file for unmatchedNormalVCF file $idx does not exist.\n") if(! -e $idx);
          $umVcfTabix = Bio::DB::HTS::Tabix->new(filename => $fileName);
        }
      }else{
        if($bedloc =~ /^(http|ftp)/) {
          $umVcfTabix = Bio::DB::HTS::Tabix->new(filename => $bedloc);
        }else{
          my $idx = $bedloc.".tbi";
          croak ("Tabix file for unmatchedNormal bed file $idx does not exist.\n") if(! -e $idx);
          $umVcfTabix = Bio::DB::HTS::Tabix->new(filename => $bedloc);
        }
      }
      $fileList->{$chr} = $umVcfTabix if(defined $umVcfTabix);
    }
  close($REF);
  return $fileList;
}

sub buildUnmatchedVCFFileListFromReference{
	my ($umLoc,$refFai,$cfg,$flag) = @_;
	my $fileList;
	my $bedloc = undef;
  #Check for bed file location before building the um normal panel
  if(defined($flag->intersect_file)){
    $bedloc = $flag->intersect_file;
  }
	if(check_exists_remote($bedloc)){
	  #Single bed file found.
    $umformat = $UNMATCHED_FORMAT_BED;
    $fileList = buildUnmatchedVCFForSingleBed($bedloc,$cfg,$refFai);
	}#End of check for unmatched bed
	else{
	  $umformat = $UNMATCHED_FORMAT_VCF;
	  $fileList = buildUnmatchedVCFForMultiVCF($bedloc,$umLoc,$cfg,$refFai,$umformat);
  }
	return $fileList;
}

sub getIntersectMatches{
	my ($vcfFile,$bedFile,$flagName) = @_;
	return 0 if(!defined $bedFile);
	die "ERROR: Unable to find $bedFile" unless(-e $bedFile);
	if(-s _ == 0) {
	  warn "WARN: Empty bed file $bedFile ... skipping intersect";
	  return 0;
	}
	#Run intersect and parse output.
	#Build command
	my $cmd = 'bedtools';
	$cmd .= ' intersect -sorted -a '.$vcfFile.' -b '.$bedFile.'';
	#Run intersect
	my $pid = open( my $IN,$cmd.' 2>&1 |') or croak("Error reading command $cmd: $?");
	while(<$IN>) {
		my $line = $_;
		chomp($line);
		if($line =~ m/ERROR/){
			croak("Problem encountered trying to intersect with $bedFile: \n",$line);
		}
		my ($chr,$pos,@unreqd) = split(/\t/,$line);
		$intersectFlagStore->{$chr.':'.$pos}->{$flagName} = 1;
	}
	close($IN) or $? != 0 or croak("Error in reading file $? (possibly trying to close reading of command $cmd. $!)");
	return;
}

sub option_builder {
	my ($factory) = @_;

	my %opts = ();

	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'i|input=s' => \$opts{'f'},
		'o|outFile=s' => \$opts{'o'},
		'c|flagConfig=s' => \$opts{'c'},
		's|species=s' => \$opts{'s'},
		't|studyType=s' => \$opts{'t'},
		'm|tumBam=s' => \$opts{'m'},
		'n|normBam=s' => \$opts{'n'},
		'g|indelBed=s' => \$opts{'g'},
		'v|flagToVcfConfig=s' => \$opts{'v'},
		'umv|unmatchedVCFLoc=s' => \$opts{'umv'},
		'ref|reference=s' => \$opts{'ref'},
		'verbose' => \$opts{'loud'},
		'l|linesToCache=s' => \$opts{'l'},
		'ab|annoBedLoc=s' => \$opts{'ab'},
		'p|processid=s' => \$opts{'p'},
		'sp|sampleToIgnoreInUnmatched=s' => \$opts{'sp'},
		'mod|module=s' => \$opts{'mod'},
		'b|bedFileLoc=s' => \$opts{'b'},
		'version' => \$opts{'version'},
	);
	return \%opts;
}

sub validateInput {
  my $opts = shift;
  pod2usage(0) if($opts->{'h'});

  if(defined $opts->{'version'}) {
    print sprintf "VERSION: %s\n", $VERSION;
    exit 1;
  }
  delete $opts->{'version'}; # needs to be deleted or breaks tests

  die( "Unknown parameter: ".$ARGV[0]) if(scalar(@ARGV) > 0);
  pod2usage("Missing parameters") if(!defined($opts->{'f'}) || !defined($opts->{'o'})
  										|| !defined($opts->{'s'})
  										|| !defined($opts->{'m'}) || !defined($opts->{'n'})
  										|| !defined($opts->{'ref'}) );


  unless(-e $opts->{'f'} && -r $opts->{'f'}){
  	pod2usage("VCF file to flag does not exist or has incorrect permissions: ".$opts->{'f'}."\n");
  }
  unless(-e $opts->{'ref'} && -r $opts->{'ref'} && $opts->{'ref'} =~ m/\.fai$/){
  	pod2usage("Reference index (fai) file doesn't exist or has bad permissions: ".$opts->{'ref'}."\n");
  }
  unless(-e $opts->{'n'} && -r $opts->{'n'}){
  	pod2usage("Non existant or unreadable normal sample bam file: ".$opts->{'n'}."\n");
  }
  unless(-e $opts->{'m'} && -r $opts->{'m'}){
  	pod2usage("Non existant or unreadable tumour sample bam file: ".$opts->{'m'}."\n");
  }

  get_config_files($opts);

  if($opts->{'g'}){
  	if(!-e $opts->{'g'} || !-r $opts->{'g'}){
  		pod2usage("Non existant or unreadable germline indel bed file: ".$opts->{'g'}."\n");
  	}
  }
  if($opts->{'umv'} && $opts->{'umv'} !~ m/^(http|ftp)/){
  	if(! -e $opts->{'umv'} || ! -d $opts->{'umv'} || ! -r $opts->{'umv'}){
  		pod2usage("Unmatched VCF location directory has bad permissions or doesn't exist: ".$opts->{'umv'}."\n");
  	}
  }
  $opts->{'s'} = uc($opts->{'s'});
  if(defined($opts->{'t'})){
    pod2usage("Unrecognised study type parameter ".$opts->{'t'}." known types: (".
                  join('|',@PERMITTED_SEQ_TYPES,).")") if(!grep(/$opts->{'t'}/i, @PERMITTED_SEQ_TYPES));
  }else{
    $opts->{'t'} = $WGS_TYPE;
    print "Using default study type of genome as no option 't' passed\n" if($opts->{'loud'});
  }
  my $type = lc($opts->{'t'});
  if(grep {/$type/} @PERMITTED_WXS_TYPES){
  	$opts->{'t'} = $WXS_TYPE;#'PULLDOWN';
  }elsif(grep {/$type/} @PERMITTED_WGS_TYPES){
  	$opts->{'t'} = $WGS_TYPE;#'GENOME';
  }elsif(grep {/$type/} @PERMITTED_RNA_TYPES){
  	$opts->{'t'} = $RNASEQ_TYPE;#'RNASEQ';
  }elsif(grep {/$type/} @PERMITTED_AMPLICON_TYPES){
  	$opts->{'t'} = $AMPLICON_TYPE;#'TARGETED';
  }else{
    croak($opts->{'t'}.": unrecognised sequencing type.");
  }
  $opts->{'t'} = uc($opts->{'t'});
  if(!defined($opts->{'l'}) || $opts->{'l'} !~ m/^\d+$/g){
  	$opts->{'l'} = $DEFAULT_LINE_CACHE;
  }
  if(defined($opts->{'b'}) && (! -e $opts->{'b'} || ! -d $opts->{'b'})){
  	pod2usage("Genome related bed file directory ".$opts->{'b'}." does not exist or is not a directory");
  }

  if(defined($opts->{'ab'}) && (! -e $opts->{'ab'} || ! -d $opts->{'ab'})){
  	pod2usage("Annotations bed file directory ".$opts->{'ab'}." does not exist or is not a directory");
  }

  if(defined($opts->{'mod'})){
    try{
      require "".$opts->{'mod'}."";
    }catch{
      pod2usage("Module provided ".$opts->{'mod'}." couldn't be used: $_");
    };
  }else{
    $opts->{'mod'} = $DEFAULT_FLAG_MODULE;
  }
  return;
}

sub get_config_files {
  my $options = shift;

  my $data_path = $Bin.'/../config/';

  $data_path = dist_dir('cgpCaVEManPostProcessing') unless(-e $data_path);

  if($options->{'v'} && (! -e $options->{'v'} || ! -r $options->{'v'})){
  	pod2usage("Error with flagToVcfConfig input check permissions.".$options->{'v'}."\n");
  }
  elsif(!$options->{'v'}){
  	$options->{'v'} = sprintf $FLAG_TO_VCF_CONFIG, $data_path;
  	print "Defaulting to use $options->{v} as flag to vcf config file.\n" if($options->{'loud'});
  	pod2usage("Default flag to vcf config file $options->{v} not found.") unless(-e $options->{'v'} && -r $options->{'v'});
  }

  if($options->{'c'} && (! -e $options->{'c'} && ! -r $options->{'c'})){
    pod2usage("Flag config file does not exist or has incorrect permissions: ".$options->{'c'}."\n");
  }
  elsif(!$options->{'c'}){
    $options->{'c'} = sprintf $FLAG_CONFIG, $data_path;
    print "Defaulting to use $options->{c} as config file.\n" if($options->{'loud'});
  	pod2usage("Default config file $options->{c} not found.") unless(-e $options->{'c'} && -r $options->{'c'});
  }

  return 1;
}

__END__

=head1 NAME

flagCaVEMan.pl - Script to flag CaVEMan VCF output.

=head1 SYNOPSIS

flagCaVEMan.pl [-h] -f vcfToFlag.vcf -o flaggedVCF.vcf -c configFile.ini -s human -t pulldown -v vcfFlagNames.ini -n norm.bam -m tum.bam [-u unmatchedStore.tmp]

  General Options:

    --help                 (-h)       Brief documentation

    --version              (-v)       Output the version number and exit

    --input                (-i)       The VCF input file to flag

    --outFile              (-o)       The VCF output file to write

    --species              (-s)       Species associated with this vcf file to use

    --tumBam               (-m)       Tumour hts (bam/cram file

    --normBam              (-n)       Normal hts (bam/cram) file

    --bedFileLoc           (-b)       Path to a folder containing the required sorted bed files (if required) i.e. the non annotation bed files.
                                      Names of files will be taken from the flag_to_vcf config file

    --indelBed             (-g)       A bed file containing germline indels to filter on

    --unmatchedVCFLoc      (-umv)     Path to a directory containing the unmatched VCF normal files listed in the
                                      config file or unmatchedNormal.bed.gz (bed file is used in preference)

    --annoBedLoc           (-ab)      Path to bed files containing annotatable regions and coding regions

    --reference            (-ref)     Reference index (fai) file corresponding to the mapping of the data being processed.

    --module               (-mod)     Path to module extending Sanger::CGP::CavemanPostProcessor::PostProcessor containing user defined flags

    --verbose

  OPTIONAL:

    --sampleToIgnoreInUnmatched    (-sp) Unmatched normal to ignore (to be used if the sample is one of those with a normal in the panel)

    --processid                    (-p)  Id anaylsis process to be added at a CGP specific header

    --flagConfig                   (-c)  Config ini file to use for flag list and settings

    --flagToVcfConfig              (-v)  Config::Inifiles style config file containing VCF flag code to flag name conversions see
                                          ../config/flag.to.vcf.convert.ini for example

    --studyType                    (-t)   Study type, used to decide parameters in file

  Examples:

    flagCaVEMan.pl [-h] -f vcfToFlag.vcf -o flaggedVCF.vcf -c configFile.ini -s human -t pulldown


=cut
