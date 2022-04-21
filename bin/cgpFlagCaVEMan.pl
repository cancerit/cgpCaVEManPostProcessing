#!/usr/bin/perl
# Copyright (c) 2014-2022
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of cgpCaVEManPostProcessing.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
#

use strict;
use warnings FATAL => 'all';

use FindBin;
use lib "$FindBin::Bin/../lib";

use Carp qw(croak cluck);
use English qw( -no_match_vars );

use Data::Dumper;

use Vcf;
use Getopt::Long qw(:config pass_through);
use Sanger::CGP::CavemanPostProcessor;
use Sanger::CGP::CavemanPostProcessor::ConfigParser;
use Sanger::CGP::CavemanPostProcessor::ExomePostProcessor;
use Sanger::CGP::CavemanPostProcessor::GenomePostProcessor;
use Sanger::CGP::CavemanPostProcessor::MNVPostProcessor;
use Pod::Usage;
use FindBin qw($Bin);
use Set::IntervalTree;
use List::Util qw (first);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Bio::DB::HTS::Tabix;
use File::Path qw(make_path);
use File::Basename;
use Const::Fast qw(const);
use LWP::Simple;
use IO::Zlib;

use File::ShareDir qw(dist_dir);

const my $FLAG_TO_VCF_CONFIG => '%s/flag.to.vcf.convert.ini';
const my $FLAG_CONFIG => '%s/%s/%s/flag.vcf.config.ini';
const my $DEFAULT_LINE_CACHE => 2000;
const my $OLD_CAVE_FLAG_BUG_FLAG => 'CB';
const my $OLD_CAVE_FLAG_BUG_DESC => 'Bug in older versions of CaVEMan means this position cannot be flagged';
const my $MNV_INFO_DESC => "Flags failed for %d base in MNV";
const my $MNV_INFO_ID => 'MNVFLAG_';
const my $GOOD_CAVE_VER => '0_3_1';

const my $CENTROMERIC_REPEAT_HIT_KEY => 'CENT';
const my $SIMPLE_REPEAT_HIT_KEY => 'SIMP';
const my $GERMLINE_INDEL_HIT_KEY => 'GID';
const my $SNP_HIT_KEY => 'SNP';
const my $ANNOTATABLE_HIT_KEY => 'ANN';
const my $CODING_ANNOTATION_HIT_KEY => 'COD';
const my $HIGH_SEQ_DEPTH_HIT_KEY => 'HSD';
const my $UNMATCHED_VCF_KEY => 'UMK';
const my $SMARTPHASE_SCORE_INFO => 'SPCONF';

const my $CENTROMERIC_REP_BED_PARAMETER => 'centromericRepeatBed';
const my $SIMPLE_REP_BED_PARAMETER => 'simpleRepeatBed';
const my $SNP_BED_PARAMETER => 'snpBed';
const my $ANNO_BED_PARAMETER => 'annotatableBed';
const my $CODING_BED_PARAMETER => 'codingBed';
const my $HSD_BED_PARAMETER => 'highSeqDepthBed';
const my $UMNORM_VCF_PARAMETER => 'unmatchedNormalVcf';
const my $UNMATCHED_VCF_LIST_PARAMETER => 'unmatchedNormalVCFList';

const my @INFO_FLAGS => qw/snpFlag codingFlag alignmentScoreReadLengthAdjustedFlag clippingMedianFlag alnScoreMedianFlag/;

const my $OLD_CAVE_BUG => 'CB';

const my $UNMATCHED_FORMAT_VCF => 'VCF';
const my $UNMATCHED_FORMAT_BED=> 'BED';

my $index = undef;
my $flagopts = undef;
my $intersectFlagStore;
my $unmatchedForOutput;
my $unmatchedFH;
my $const = 'Sanger::CGP::CavemanPostProcessor::Constants';
my $umformat = $UNMATCHED_FORMAT_VCF;

my $opts = option_builder();
validateInput($opts);
my $checkSize = 0;

if(defined($opts->{'idx'})){
  $index = $opts->{'idx'};
}

my @lineCache = ();

my $oldCaveVersion = 0;

if(defined($index) && $index > 0){
  $opts->{'f'} .= ".".$index;
  $opts->{'o'} .= ".".$index;
}

warn "Working on ".$opts->{'f'}." and ".$opts->{'o'}."\n" if($opts->{'loud'});

eval{
  main($opts);
};
if($EVAL_ERROR){
  croak($EVAL_ERROR);
}

#####################
#                   #
#    SUBROUTINES    #
#                   #
#####################

sub _tabix_to_interval_tree {
  my $bed = shift;
  my %tree;
  my $z = IO::Uncompress::Gunzip->new($bed, MultiStream => 1) or die "gunzip failed: $GunzipError\n";
  my $value = 1;
  while(my $line = <$z>) {
    next if ($line =~ m/^#/);
    chomp $line;
    my ($chr, $s, $e) = split /\t/, $line;
    $tree{$chr} = Set::IntervalTree->new() unless(exists $tree{$chr});
    $tree{$chr}->insert(\$value, $s, $e);
  }
  close $z;
  return \%tree;
}

sub main{
  my ($opts) = @_;
  #Get flag list and config for flagger from config.ini file
  my $umNormVcf;
  my $tabixList;
  my ($configParams,$flagList,$centBed,$simpBed,$snpBed,
        $indelBed,$annoBed,$codingBed,$hsdBed,$mnvflaglist) = setupFromConfig($opts);
  warn "Performing intersects\n" if($opts->{'loud'});
    if(grep(/centromericRepeatFlag/,@$flagList)){
      my $idx = $centBed.".tbi";
      #Run intersects for each of the potential flags, will skip if there's no file.
      croak ("Tabix file for centromericRepeatFlag file $idx does not exist.\n") if(! -e $idx);
            $tabixList->{'centromericRepeatFlag'} = _tabix_to_interval_tree($centBed);
    }
    if(grep(/simpleRepeatFlag/,@$flagList)){
      my $idx = $simpBed.".tbi";
      croak ("Tabix file for simpleRepeatFlag file $idx does not exist.\n") if(! -e $idx);
            $tabixList->{'simpleRepeatFlag'} = _tabix_to_interval_tree($simpBed);
    }
    if(grep(/snpFlag/,@$flagList)){
      my $idx = $snpBed.".tbi";
      croak ("Tabix file for snpFlag file $idx does not exist.\n") if(! -e $idx);
      $tabixList->{'snpFlag'} = Bio::DB::HTS::Tabix->new(filename => $snpBed);
    }
    if(grep(/germlineIndelFlag/,@$flagList)){
      my $idx = $indelBed.".tbi";
      croak ("Tabix file for germlineIndelFlag file $idx does not exist.\n") if(! -e $idx);
      $tabixList->{'germlineIndelFlag'} = _tabix_to_interval_tree($indelBed);
    }
    if(grep(/annotationFlag/,@$flagList)){
      my $idx = $annoBed.".tbi";
      croak ("Tabix file for annotationFlag file $idx does not exist.\n") if(! -e $idx);
      $tabixList->{'annotationFlag'} = _tabix_to_interval_tree($annoBed);
    }
    if(grep(/codingFlag/,@$flagList)){
      my $idx = $codingBed.".tbi";
      croak ("Tabix file for codingFlag file $idx does not exist.\n") if(! -e $idx);
      $tabixList->{'codingFlag'} = _tabix_to_interval_tree($codingBed);
    }
    if(grep(/hiSeqDepthFlag/,@$flagList)){
      warn "Performing HSD intersect\n" if($opts->{'loud'});
      my $idx = $hsdBed.".tbi";
      croak ("Tabix file for hiSeqDepthFlag file $idx does not exist.\n") if(! -e $idx);
      $tabixList->{'hiSeqDepthFlag'} = _tabix_to_interval_tree($hsdBed);
    }
    #Setup postprocessing module
    my $unmatchedVCFFlag = 0;
    if(grep(/unmatchedNormalVcfFlag/,@$flagList)){
      warn "Generating hash of unmatched normal VCF files\n" if($opts->{'loud'});
      $unmatchedVCFFlag = 1;
      #Build up a list of VCF files present in the location passed, there will be on per chromosome....
      $umNormVcf = buildUnmatchedVCFFileListFromReference($opts->{'umv'},$opts->{'ref'},$configParams);
    }
  #Setup postprocessing module
  #my $flagger;
  my $flagger = initFlagModuleForSpeciesType($configParams,$opts->{'t'});
  my $VCFOUT;
  my $infoWritten = 0;
  #Open VCF input file
  my $vcf = Vcf->new(file=>$opts->{'f'}, version=>'4.1');
  warn "Starting flagging\n" if($opts->{'loud'});

  # Get full path of output file.
  my ($fname,$dirs,$suffix) = fileparse($opts->{'o'});
  # If the directory doesn't exist. Create it.
  if(! -e $dirs){
    make_path($dirs) or croak("Error trying to create output directory '$dirs'");
  }

  # Open output VCF file
  open($VCFOUT , '>', $opts->{'o'}) || croak("Error trying to open VCF output file ".$opts->{'o'}.": $!");
    $vcf->parse_header();
    #$vcf->format_header();
    #Get the CaVEMan version from the header so we know if we need the SWBad flag included.
    my $caveVer = getCavemanVersionFromVCF($vcf);
    $oldCaveVersion = isOldVersionOfCaVEMan($caveVer);
    #load flagtoVCF config ini
    my $cfg = Config::IniFiles->new( -file => $opts->{'v'} );

    if(exists($opts->{'p'}) && defined($opts->{'p'})){
      $vcf->add_header_line({key=>'cgpAnalysisProc',value=>$opts->{'p'}}, append=>1);
    }
    #add flag list for this species/study type combination to FILTER fields.
    appendFiltersToHeader($vcf,$flagList,$cfg,$configParams,$oldCaveVersion);
    #If our flaglist contains any of the bedtools intersect required flags we nbeed those to run over the whole file.
    #runBedtoolsForFlags($opts->{'f'},$flagList);
    foreach my $key(keys %$opts){
      if($key =~ m/^[focmngvub]|(db|um|ab|ref|umv)$/ && defined($opts->{$key})){
        $opts->{$key} = _trim_file_path($opts->{$key});
      }
      if(!defined($opts->{$key})){
        $opts->{$key} = '.';
      }
    }
    $vcf->add_header_line({
              key=>'vcfProcessLog',
              'InputVCF'=>''._trim_file_path($opts->{'f'}).'',
              'InputVCFSource'=>'FlagCaVEManVCF.pl',
              'InputVCFVer'=>''.get_version().'',
              'InputVCFParam'=>$opts
              } ,append => 1);

    print $VCFOUT $vcf->format_header();
    #Iterate through each VCF variant line
    while (my $x=$vcf->next_data_array()){
      my $isInUmVCF = undef;
      #Clear the flags and SNP status before flagging.
      $$x[6]='.'; #This resets all flags.
      #This ensures we haven't got any info flags (SNP) in existance...
      $$x[7]=$vcf->add_info_field($$x[7],'SNP'=>undef,'coding'=>undef,'ASRD'=>undef,
                                                        'CLPM'=>undef,'ASMD'=>undef,'SR'=>undef);
      #Files are sorted by chr/pos, so we can open the tabix index once per vcf file.
      my $this_flagger;
      my $results = ();
      my $results_arr = [];
      if (length($$x[3]) > 1) {
        $this_flagger = $flagger->{'MNV'};
        #Divide MNV into SNV sections and flag each as an individual SNV
        my $mnv_len = length($$x[3])-1;
        for my $i(0 .. $mnv_len){
          my $thisPos = $$x[1] + $i;
          if(defined($mnvflaglist) && first {'unmatchedNormalVcfFlag'} @$mnvflaglist && exists($umNormVcf->{$$x[0]})){
            $isInUmVCF = getUnmatchedVCFIntersectMatch($$x[0],$thisPos,$umNormVcf->{$$x[0]},$UNMATCHED_VCF_KEY);
          }
          #Use $mnvflaglist of MNV permitted flags instead of $flagList
          $results = getVCFToAddResultsOfFilters($$x[0],$thisPos,substr($$x[3],$i,1),substr($$x[4],$i,1),$mnvflaglist,$this_flagger,$cfg,$x,$vcf,$configParams,$isInUmVCF,$tabixList,$i);
          push(@$results_arr,$results);
        }
      }else{
        if($unmatchedVCFFlag==1 && exists($umNormVcf->{$$x[0]})){
          $isInUmVCF = getUnmatchedVCFIntersectMatch($$x[0],$$x[1],$umNormVcf->{$$x[0]},$UNMATCHED_VCF_KEY);
        }
        $this_flagger = $flagger->{'SNV'};
        $results = getVCFToAddResultsOfFilters($$x[0],$$x[1],$$x[3],$$x[4],$flagList,$this_flagger,$cfg,$x,$vcf,$configParams,$isInUmVCF,$tabixList,$index);
      }
      #Only if we have a list of mnv usable flags and the MNVflag itself
      my $mnv_flagname = $const->flag_names('mnvFlag');
      if(defined $mnvflaglist && scalar(@$mnvflaglist) > 0 && first {$mnv_flagname} @$flagList && length($$x[3]) > 1){
        #TODO run mnv flag here, use list of mnv flags and get results of all flags to apply to each position in the MNV
        #Run MNV specific flag here
        ($results, $vcf, $x) = $flagger->{'MNV'}->run_mnv_flag($results_arr,$vcf,$x);
      }
      #Add the relevant filters or PASS to the filter section.
      $$x[6]=$vcf->add_filter($$x[6],%$results);
      #Validate this line and append this line to the output file;
      push(@lineCache,correct_flag_sort($vcf->format_line($x)));
      if(scalar(@lineCache)>=$opts->{'l'}){
        foreach my $line(@lineCache){
          print $VCFOUT $line;
        }
        @lineCache = ();
      }
    }#End of iteration throuh each position
    #Lastly empty the cache of vcf lines to output.
    foreach my $line(@lineCache){
      print $VCFOUT $line;
    }
    @lineCache = ();
  close($VCFOUT);

  warn "Done flagging\n" if($opts->{'loud'});
}

sub close_tabix{
  my ($tabixList) = @_;
  foreach my $flagType(keys %$tabixList){
    if(defined($tabixList->{$flagType})){
      $tabixList->{$flagType}->close;
    }
  }
  return;
}

sub correct_flag_sort{
  my ($line) = @_;
  my @all_fields = split /\t/, $line;
  my @flags = split /;/,$all_fields[6];
  $all_fields[6] = join ";" , sort {$a cmp $b} @flags;
  $line = join "\t", @all_fields;
  return $line;
}

sub check_exists_remote{
  my ($fileloc) = @_;
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

sub _checkConfigAndBedUnmatched{
  my ($bedloc,$cfg) = @_;
  my $FH;
  my $check = 0;
  $FH = new IO::Zlib;
  $FH->open($bedloc,"rb") or croak("Error opening zipped bed file '$bedloc': $!");
    while(<$FH>){
      my $line = $_;
      if($line =~ m/^\s*#PARAMS/){
        #Unmatched normal panel for CaVEMan
        #PARAMS vcfUnmatchedMinMutAlleleCvg=3   vcfUnmatchedMinSamplePct=1.000  panelMD5=f7150610e2af583780f55cc58e278063
        chomp($line);
        my ($tmp,$tmp2,$vcfminmut,$tmp3,$vcfminsamp,$tmp4,$md5) = split /\s+|=/, $line ;
        if($vcfminmut == $cfg->{"vcfUnmatchedMinMutAlleleCvg"} && $vcfminsamp ==  $cfg->{"vcfUnmatchedMinSamplePct"}){
          $check = 1;
          last;
        }
      }elsif($line !~ m/^\s*#/){#Assumes header line comes first. This will cause an error if it's not or it's not present.
        last;
      }
    }
  $FH->close or croak("Error closing zipped bed file '$bedloc': $!");
  return $check;
}

sub buildUnmatchedVCFFileListFromReference{
  my ($umLoc,$refFai,$cfg) = @_;
  my $fileList;
  my $REF;
  my $count = 0;
  my $bedloc = $umLoc."/unmatchedNormal.bed.gz";
  if(check_exists_remote($bedloc)){
    #Single bed file found.
    $umformat = $UNMATCHED_FORMAT_BED;
    #Check here that the header matches the params set in the configs.
    croak("Unmatched panel parameters in bed provided don't match those in the config file.") if(_checkConfigAndBedUnmatched($bedloc,$cfg)!=1);
    #Check the tabix file exists.
    my $umBedTabix = $bedloc.".tbi";
    croak("Unmatched bedfile $umBedTabix") if (! -e $umBedTabix);
    my $umBedTbi = Bio::DB::HTS::Tabix->new(filename => $bedloc);
    open($REF, '<', $refFai) or croak("Error opening reference index $refFai: $!");
      while(<$REF>){
        my $line = $_;
        next if($line =~ m/^\s*#/);
        chomp($line);
        my ($chr,undef) = split(/\t/,$line);
        $fileList->{$chr} = $umBedTbi;
      }
    close($REF);
  }
  else{
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
  }
  return $fileList;
}

sub _interval_hit {
    my ($tabix, $chr, $from, $to) = @_;
    return 0 unless (exists $tabix->{$chr});
    return scalar @{$tabix->{$chr}->fetch($from-1, $to)};
}

sub _interval_hit_data {
    my ($tabix, $chr, $from, $to) = @_;
    return 0 unless (exists $tabix->{$chr});
    return $tabix->{$chr}->fetch($from-1, $to);
}

sub isOldVersionOfCaVEMan{
  my ($ver) = @_;
  return 1 if(!defined($ver) || $ver eq "");
  my ($good_main,$good_feat,$good_bug) = split(/[_\.]/,$GOOD_CAVE_VER);
  $ver =~ s/^[^\|]+\|//g;
  $ver =~ s/\|[^\|]+$//g;
  my ($main,$feat,$bug) = split(/[\._]/,$ver);
  return 1 if( $main < $good_main || ($main == $good_main && $feat < $good_feat) || ($main == $good_main && $feat == $good_feat && $bug < $good_bug));
  return 0;
}

sub getCavemanVersionFromVCF{
  my ($vcf) = @_;
  my $caveVer = undef;
  my $verLine = $vcf->get_header_line(key=>'cavemanVersion');
  $caveVer = $verLine->[0]->[0]->{'value'};
  #This might be an original CaVEMan VCF or an old run before the headers were fixed
  if(!defined($caveVer) || $caveVer eq ""){
    my $plLines = $vcf->get_header_line(key=>'vcfProcessLog');
    foreach my $plln(@$plLines){
      if($plln->[0]->{'InputVCFSource'} eq 'CaVEMan'){
        $caveVer = $plln->[0]->{'InputVCFVer'};
      }
    }
  }
  #Lastly, perhaps the VCF version is old so we need to cycle through each header ourselves to find the right one :-/.
  if(!defined($caveVer) || $caveVer eq ""){
    my $head = $vcf->format_header();
    my @plLines = split(/\n/,$head);
    foreach my $plln(@plLines){
      next unless($plln =~ m/^\s*#+vcfProcessLog.+InputVCFSource=<CaVEMan>.+/);
      $plln =~ s/^\s*#+//g;
      my @ln = split(/=/,$plln,2);
      my $ky = $ln[0];
      my @rest = split(/>,/,$ln[1]);
      #Look through the body to find the version
      foreach my $rt(@rest){
        next unless $rt =~ m/InputVCFVer/;
        $rt =~ s/<//g;
        my $tmp;
        (undef,$caveVer) = split(/=/,$rt);
        last;
      }
    }
  }
  if(!defined($caveVer) || $caveVer eq ""){
    $caveVer = ".";
  }  return $caveVer;
}

sub _trim_file_path{
  my ($string) = @_;
  my @bits = (split("/", $string));
  return pop @bits;
}

sub getVCFToAddResultsOfFilters{
  my ($chr,$pos,$wt,$mut,$flagList,$flagger,$cfg,$x,$vcf,$configParams,$isInUmVCF,$tabixList,$index) = @_;
  my %resFlags = ();
  $flagger->runProcess($chr,$pos,$pos,$wt,$mut);
  if($oldCaveVersion && (!exists($flagger->_muts->{'tqs'}) || scalar(@{$flagger->_muts->{'tqs'}})== 0)){
    $resFlags{$OLD_CAVE_BUG} = 1;
    return \%resFlags;
  }
  #Iterate through each flag we need to run.
  foreach my $flagName(@$flagList){
    next if $flagName eq 'mnvFlag';
    #Get the id of this flag.
    my $flagId = $cfg->val($flagName,"id");
    #Run this flag
    my $flagRes = runFlagger($flagger,$flagName,$flagId,$chr,$pos,$mut,$x,$vcf,$configParams,$isInUmVCF,$tabixList,$index);
    #Add inverse of the flag result to the resFlags (we return 1 for pass... fails want to be added to filter)
    if($flagName !~ m/(snp|coding)Flag/){
      $resFlags{$flagId} = !$flagRes;
    }
  }
  return \%resFlags;
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
  my $IN;
  my $pid = open($IN,$cmd.' 2>&1 |');
  while(<$IN>) {
    my $line = $_;
    chomp($line);
    if($line =~ m/ERROR/){
      croak("Problem encountered trying to intersect with $bedFile: \n",$line);
    }
    my ($chr,$pos,@unreqd) = split(/\t/,$line);
    $intersectFlagStore->{$chr.':'.$pos}->{$flagName} = 1;
  }
  close($IN);
  return;
}

sub getUnmatchedVCFIntersectMatch{
  my ($chr,$pos,$tabix,$flagName) = @_;
  #Only check for positions we've not already sorted.
  # new tabix is 1-based for both coordinates
  my $iter = $tabix->query_full($chr,$pos,$pos);
  my $line = undef;
  $line = $iter->next if(defined($iter)); # undef if not found
  return $line;
}

sub runFlagger{
  my ($flagger,$flagName,$flagId,$chr,$pos,$mut,$x,$vcf,$cfg,$isInUmVCF,$tabixList,$index) = @_;
  my $coord = $chr.':'.$pos;
  if($flagName eq 'depthFlag'){
    #DEPTH
    return $flagger->getDepthResult();
  }elsif($flagName eq 'readPositionFlag'){
    #READ POSITION
    return $flagger->getReadPositionResult();
  }elsif($flagName eq 'matchedNormalFlag'){
    #MATCHED NORMAL
    return $flagger->getNormMutsAllelesResult();
  }elsif($flagName eq 'pentamericMotifFlag'){
    #PENTAMERIC MOTIF
    return $flagger->getPentamerResult();
  }elsif($flagName eq 'avgMapQualFlag'){
    #MEAN MAP QUAL
    return $flagger->getAvgMapQualResult();
  }elsif($flagName eq 'withinGapRangeFlag'){
    return $flagger->getReadGapFlagResult();
  }elsif($flagName eq 'germlineIndelFlag'){
    #GERMLINE INDEL
    return !_interval_hit($tabixList->{$flagName},$chr,$pos,$pos);
  }elsif($flagName eq 'tumIndelDepthFlag'){
    #TUM INDEL DEPTH
    return $flagger->getTumIndelReadDepthResult();
  }elsif($flagName eq 'sameReadPosFlag'){
    #RD POS STACKING
    return $flagger->getDifferingReadPositionResult();
  }elsif($flagName eq 'simpleRepeatFlag'){
    #SIMPLE REPEATS
        return !_interval_hit($tabixList->{$flagName},$chr,$pos,$pos);
  }elsif($flagName eq 'unmatchedNormalVcfFlag'){
    if(defined($isInUmVCF)){
      if($umformat eq $UNMATCHED_FORMAT_BED){#Check for bed rather than VCF

            my ($chr,$start,$stop,$score,$mutlist) = split(/\t/,$isInUmVCF);
            return 0 if($mutlist =~ m/$mut/i);

            #End of if this is BED unmatched
      }elsif($umformat eq $UNMATCHED_FORMAT_VCF){

        my ($ch,$po,$ident,$refAll,$altAll,$quality,$filts,$info,$format,@samples) = split(/\t/,$isInUmVCF);
        #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE1b  SAMPLE2b  SAMPLE3b
        my ($geno,@formats) = split(':',$format);
        my $sampleHitCount = 0;
        my $totalSampleCnt = 0;
        foreach my $sampData(@samples){
          $totalSampleCnt++;
          next if($sampData eq q{0:.:.:.:.} || $sampData eq q{-} || $sampData eq q{.});
          #GT:GF:CF:TF:AF  0|0:41:0:0:0
          my ($gentype,@data) = split(':',$sampData);
          my $totalCvg = 0;
          my $mutAlleleCvg = 0;
          my $proportion = 0;
          for (my $i=0;$i<scalar(@data);$i++){
            next unless ($formats[$i] =~ m/^[ACGT]{1}[A-Z]$/);
            $totalCvg += $data[$i];
            if(substr($formats[$i],0,1) eq $mut){
              $mutAlleleCvg = $data[$i];
            }
          }
          if($mutAlleleCvg >= $cfg->{"vcfUnmatchedMinMutAlleleCvg"}){
            $sampleHitCount++;
          }
        }
        return 1 if($totalSampleCnt == 0);
        if((($sampleHitCount/$totalSampleCnt)*100) >= $cfg->{"vcfUnmatchedMinSamplePct"}){
          return 0;
        }
      }#End of if this is VCF

    }
    return 1;
  }elsif($flagName eq 'centromericRepeatFlag'){
    #CENTROMERIC REPEATS
    #Use intersect to check for centromeric repeats
    return !_interval_hit($tabixList->{$flagName},$chr,$pos,$pos);
  }elsif($flagName eq 'snpFlag'){
    #SNPS
    my $iter = $tabixList->{$flagName}->query_full($chr,$pos,$pos);
    my $line = undef;
    $line = $iter->next if(defined($iter)); # undef if not found
    if(defined($line)){
      $$x[7]=$vcf->add_info_field($$x[7],$flagId=>'');
    }
    return -1;
  }elsif($flagName eq 'cavemanMatchNormalProportionFlag'){
    return $flagger->getCavemanMatchedNormalResult($vcf,$x,$index);
  }elsif($flagName eq 'phasingFlag'){
    #PHASING
    return $flagger->getPhasingResult();
  }elsif($flagName eq 'annotationFlag'){
    #ANNOTATION
    return _interval_hit($tabixList->{$flagName},$chr,$pos,$pos);
  }elsif($flagName eq 'hiSeqDepthFlag'){
    #HIGH SEQ DEPTH
    return !_interval_hit($tabixList->{$flagName},$chr,$pos,$pos);
  }elsif($flagName eq 'codingFlag'){
    #CODING
    if(_interval_hit($tabixList->{$flagName},$chr,$pos,$pos)){
      $$x[7]=$vcf->add_info_field($$x[7],$flagId=>'');
    }
    return -1;
  }elsif($flagName eq 'clippingMedianFlag'){
    $$x[7]=$vcf->add_info_field($$x[7],$flagId=>$flagger->getClipMedianResult());
  }elsif($flagName eq 'alignmentScoreReadLengthAdjustedFlag'){
    $$x[7]=$vcf->add_info_field($$x[7],$flagId=>$flagger->getAlignmentScoreMedianReadAdjusted());
  }elsif($flagName eq 'alnScoreMedianFlag'){
    $$x[7]=$vcf->add_info_field($$x[7],$flagId=>$flagger->getAlignmentScoreMedian());
  }elsif($flagName eq 'lowMutBurdenFlag'){
    #Low mut burden
    croak('LOW_MUT_BURDEN');
  }elsif($flagName eq 'singleEndFlag'){
    #SINGLE STRAND
    return $flagger->getSingleEndResult();
  }elsif($flagName eq 'matchedNormalProportion'){
    return $flagger->getMatchedNormalProportionResult();
  }else{
    croak("Unrecognised flag name $flagName\n");
  }
}

sub appendFiltersToHeader{
  my ($vcf,$flagList,$cfg,$configParams,$oldCave) = @_;
  foreach my $flagName(@$flagList){
    if(!$cfg->val($flagName,"id") || !$cfg->val($flagName,"description")){
      croak("No config description or ID found for flag ".$flagName."\n");
    }
    my $desc = $cfg->val($flagName,"description");
    my $info = $cfg->val($flagName,"info");
    my $newDesc = getDescriptionWithParams($desc,$configParams);
    if($info==0){
      $vcf->add_header_line({key=>'FILTER', ID=>$cfg->val($flagName,"id"),Description=>$newDesc});
    }else{
      my $val = $cfg->val($flagName,"val");
      my $type = $cfg->val($flagName,"type");
      $vcf->add_header_line({key=>'INFO', ID=>$cfg->val($flagName,"id"),Type=>$type,Number=>$val,Description=>$newDesc});
    }
  }
  if($oldCave){
    $vcf->add_header_line({key=>'FILTER', ID=>$OLD_CAVE_FLAG_BUG_FLAG,Description=>$OLD_CAVE_FLAG_BUG_DESC});
  }
  #See if we have (if any) FORMAT lines for MNVs to find the max length
  my $max_len = 1;
  while (defined($vcf->get_header_line(key=>'FORMAT', ID=>'FAZ_'.($max_len+1))) &&
                scalar(@{$vcf->get_header_line(key=>'FORMAT', ID=>'FAZ_'.($max_len+1))})){
    $max_len += 1;
  }
  if ($max_len > 1){ # We have MNVs
    my $type = "String";
    my $val = ".";
    for my $i(1..$max_len){
      my $description = sprintf($MNV_INFO_DESC,$i);
      my $id = "$MNV_INFO_ID"."$i";
      $vcf->add_header_line({key=>'INFO', ID=>$id,Type=>$type,Number=>$val,Description=>$description});
    }
  }
  return;
}

sub getDescriptionWithParams{
  my ($desc,$configParams) = @_;
  foreach my $pName(keys %$configParams){
    my $val = $configParams->{$pName};
    $desc =~ s/$pName/$val/g;
  }
  return $desc;
}

sub setupFromConfig{
  my ($opts) = @_;
  #Load in config file.
  my ($configParams,$flagList,$bedFileParams, $mnvflaglist) = Sanger::CGP::CavemanPostProcessor::ConfigParser::getConfigParams($opts,$flagopts);
  #Use the parameters and do some sanity checks for flag/bed file requirements.
  my ($centBed,$simpBed,$snpBed,$indelBed,$annoBed,$codingBed,$hsdBed) = undef;
  my @errs = ();
  if(grep(/centromericRepeatFlag/,@$flagList)){
    if(!defined($bedFileParams->{$CENTROMERIC_REP_BED_PARAMETER})){
      push(@errs,"Found no file path to centromeric repeat bed file. Set ".$CENTROMERIC_REP_BED_PARAMETER." in the config input\n");
    }
    $centBed = $opts->{'b'}.'/'.$bedFileParams->{$CENTROMERIC_REP_BED_PARAMETER};
    if(! -e $centBed){
      push(@errs,"Specified bed file does not exist: $centBed\n");
    }
  }

  if(grep(/simpleRepeatFlag/,@$flagList)){
    if(!defined($bedFileParams->{$SIMPLE_REP_BED_PARAMETER})){
      push(@errs,"Found no file path to simple repeat bed file. Set ".$SIMPLE_REP_BED_PARAMETER." in the config input\n");
    }
    $simpBed = $opts->{'b'}.'/'.$bedFileParams->{$SIMPLE_REP_BED_PARAMETER};
    if(! -e $simpBed){
      push(@errs,"Specified bed file does not exist: $simpBed\n");
    }
  }

  if(grep(/snpFlag/,@$flagList)){
    if(!defined($bedFileParams->{$SNP_BED_PARAMETER})){
      push(@errs,"Found no file path to SNP bed file. Set ".$SNP_BED_PARAMETER." in the config input\n");
    }
    $snpBed = $opts->{'b'}.'/'.$bedFileParams->{$SNP_BED_PARAMETER};
    if(! -e $snpBed){
      push(@errs,"Specified bed file does not exist: $snpBed\n");
    }
  }

  if(grep(/germlineIndelFlag/,@$flagList)){
    if(!defined($opts->{'g'})){
      push(@errs,"Found no file path to germline indel bed file. Set -g at the commandline\n");
    }
    $indelBed = $opts->{'g'};
    if(! -e $indelBed){
      push(@errs,"Specified bed file does not exist: $indelBed\n");
    }
  }

  if(grep(/annotationFlag/,@$flagList)){
    if(!defined($bedFileParams->{$ANNO_BED_PARAMETER})){
      push(@errs,"Found no file path to annotatable bed file. Set ".$ANNO_BED_PARAMETER." in the config input\n");
    }
    $annoBed = $opts->{'ab'}.'/'.$bedFileParams->{$ANNO_BED_PARAMETER};
    if(! -e $annoBed){
      push(@errs,"Specified bed file does not exist: $annoBed\n");
    }
  }

  if(grep(/codingFlag/,@$flagList)){
    if(!defined($bedFileParams->{$CODING_BED_PARAMETER})){
      push(@errs,"Found no file path to coding annotation bed file. Set ".$CODING_BED_PARAMETER." in the config input\n");
    }
    $codingBed = $opts->{'ab'}.'/'.$bedFileParams->{$CODING_BED_PARAMETER};
    if(! -e $codingBed){
      push(@errs,"Specified bed file does not exist: $codingBed\n");
    }
  }

  if(grep(/hiSeqDepthFlag/,@$flagList)){
    if(!defined($bedFileParams->{$HSD_BED_PARAMETER})){
      push(@errs,"Found no file path to high seq depth bed file. Set ".$HSD_BED_PARAMETER." in the config input\n");
    }
    $hsdBed = $opts->{'b'}.'/'.$bedFileParams->{$HSD_BED_PARAMETER};
    if(! -e $hsdBed){
      push(@errs,"Specified bed file does not exist: $hsdBed\n");
    }
  }
  if(scalar(@errs) > 0){
    my $err = join("",@errs);
    croak("Errors encountered setting up locations from config files:\n".$err."\n");
  }
  #Add bam/cram files to config
  $configParams->{'tumBam'} = $opts->{'m'};
  $configParams->{'normBam'} = $opts->{'n'};
  #Add ref so cram works too
  my $tmp_ref = $opts->{'ref'};
  $tmp_ref =~ s/\.fai$//;
  $configParams->{'ref'} = $tmp_ref;
  return ($configParams,$flagList,$centBed,$simpBed,$snpBed,$indelBed,$annoBed,$codingBed,$hsdBed, $mnvflaglist);
}

sub initFlagModuleForSpeciesType{
  my ($params,$type) = @_;
  my $flag_mods;
  if(lc($type) eq lc("pulldown") || lc($type) eq lc("followup") || $type eq "WXS" || $type eq "AMPLICON" || $type eq "TARGETED"){
    $flag_mods->{'SNV'} = Sanger::CGP::CavemanPostProcessor::ExomePostProcessor->new(%$params);
  }elsif($type eq "WGS" || $type eq "RNASEQ"){
    $flag_mods->{'SNV'} = Sanger::CGP::CavemanPostProcessor::GenomePostProcessor->new(%$params);
  }else{
    croak("No flagging module for type: $type\n");
  }
  $flag_mods->{'MNV'} = Sanger::CGP::CavemanPostProcessor::MNVPostProcessor->new(%$params);
  return $flag_mods;
}

sub get_version {
  return Sanger::CGP::CavemanPostProcessor->VERSION;
}

sub get_config_files {
  my $options = shift;

  my $data_path = $Bin.'/../config/';

  print "No developer environment found, checking for config files passed at commandline\n" if($options->{'loud'});

  $data_path = dist_dir('cgpCaVEManPostProcessing') unless(-e $data_path);

  print "Using $data_path as share directory if files not pathed at commandline\n" if($options->{'loud'});

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
    $options->{'c'} = sprintf $FLAG_CONFIG, $data_path, $options->{'s'}, $options->{'sa'};
    print "Attempting to use $options->{c} as config file.\n" if($options->{'loud'});
    pod2usage("Default config file $options->{c} not found.") unless(-e $options->{'c'} && -r $options->{'c'});
  }

  return 1;
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
    'sa|species-assembly=s' => \$opts{'sa'},
    't|studyType=s' => \$opts{'t'},
    'm|tumBam=s' => \$opts{'m'},
    'n|normBam=s' => \$opts{'n'},
    'g|indelBed=s' => \$opts{'g'},
    'v|flagToVcfConfig=s' => \$opts{'v'},
    'umv|unmatchedVCFLoc=s' => \$opts{'umv'},
    'ref|reference=s' => \$opts{'ref'},
    'idx|index=s' => \$opts{'idx'},
    'verbose' => \$opts{'loud'},
    'l|linesToCache=s' => \$opts{'l'},
    'ab|annoBedLoc=s' => \$opts{'ab'},
    'p|processid=s' => \$opts{'p'},
    'sp|sampleToIgnoreInUnmatched=s' => \$opts{'sp'},
    'b|bedFileLoc=s' => \$opts{'b'},
    'f|flags=s' => \@{$opts{'flags'}},
    'r|flag-mnv' => \$opts{'flagmnv'},
    'version' => \$opts{'version'},
  );
  return \%opts;
}

sub validateInput {
  my $opts = shift;
  pod2usage(0) if($opts->{'h'});

  if(defined $opts->{'version'}) {
    print sprintf "VERSION: %s\n", Sanger::CGP::CavemanPostProcessor->VERSION;
    exit 0;
  }
  delete $opts->{'version'}; # needs to be deleted or breaks tests

  die( "Unknown parameter: ".$ARGV[0]) if(scalar(@ARGV) > 0);
  pod2usage("Missing parameter 'input'") if(!defined($opts->{'f'}));
  pod2usage("Missing parameter 'outFile'") if(!defined($opts->{'o'}));
  pod2usage("Missing parameter 'species'") if(!defined($opts->{'s'}));
  pod2usage("Missing parameter 'tumBam'") if(!defined($opts->{'m'}));
  pod2usage("Missing parameter 'normBam'") if(!defined($opts->{'n'}));
  pod2usage("Missing parameter 'reference'") if(!defined($opts->{'ref'}));
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

  if(defined($opts->{'flags'}) && scalar(@{$opts->{'flags'}})>0){
    @$flagopts =  split(/,/,join(',',@{$opts->{'flags'}}));
  }
  delete $opts->{'flags'};

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
    pod2usage("Unrecognised study type parameter ".$opts->{'t'}.
          " known types: (pulldown|exome|WGS|WXS|AMPLICON|genome|genomic|followup|targeted|targetted|rna_seq|rna-seq|rnaseq)") if($opts->{'t'} !~
                                m/^(pulldown|exome|WGS|WXS|AMPLICON|genome|genomic|followup|targeted|targetted|rna_seq|rna-seq|rnaseq)$/i);
  }else{
    $opts->{'t'} = 'WGS';
    print "Using default study type of WGS as not option t passed\n" if($opts->{'loud'});
  }

  if(uc($opts->{'t'}) eq 'EXOME' || uc($opts->{'t'}) eq 'PULLDOWN' || uc($opts->{'t'}) eq 'WXS'){
    $opts->{'t'} = 'WXS';
  }elsif(uc($opts->{'t'}) eq 'GENOMIC' || uc($opts->{'t'}) eq 'WGS' || uc($opts->{'t'}) eq 'GENOME'){
    $opts->{'t'} = 'WGS';
  }elsif(uc($opts->{'t'}) eq 'RNA_SEQ'||uc($opts->{'t'}) eq 'RNA-SEQ'){
    $opts->{'t'} = 'RNASEQ';
  }elsif(uc($opts->{'t'}) eq 'TARGETED'){
    $opts->{'t'} = 'TARGETED';
  }elsif(uc($opts->{'t'}) eq 'AMPLICON' || uc($opts->{'t'}) eq 'FOLLOWUP'){
    $opts->{'t'} = 'AMPLICON';
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

  pod2usage("Index is not numeric and > 0") if(defined($opts->{'idx'}) && ($opts->{'idx'} !~ m/^\d+$/g || $opts->{'idx'} < 1 ));
  return;
}

__END__

=head1 NAME

cgpFlagCaVEMan.pl - An Index aware copy of FlagCaVEManVCF.pl script and appends index to the input/output files.

=head1 SYNOPSIS

cgpFlagCaVEMan.pl [-h] -f vcfToFlag.vcf -o flaggedVCF.vcf -c configFile.ini -s human -t pulldown -v vcfFlagNames.ini -n norm.bam -m tum.bam [-u unmatchedStore.tmp]

  General Options:

    --help                 (-h)       Brief documentation

    --version              (-version) Output the version number and exit

    --input                (-i)       The VCF input file to flag.

    --outFile              (-o)       The VCF output file to write.

    --species              (-s)       Species associated with this vcf file to use.

    --species-assembly     (-sa)      Species assembly for (output in VCF)

    --tumBam               (-m)       Tumour bam file

    --normBam              (-n)       Normal bam file

    --bedFileLoc           (-b)       Path to a folder containing the centromeric, snp, hi sequence depth,
                                      and simple repeat sorted (gzipped and tabixed) bed files (if required) i.e. the non annotation bed files.
                                      Names of files will be taken from the config file.

    --indelBed             (-g)       A bed file containing germline indels to filter on

    --unmatchedVCFLoc      (-umv)     Path to a directory containing the unmatched VCF normal files listed in the
                                      config file or unmatchedNormal.bed.gz (bed file is used in preference).

    --annoBedLoc           (-ab)      Path to bed files containing annotatable regions and coding regions.

    --reference            (-ref)     Reference index (fai) file corresponding to the mapping of the data being processed.
                                        (must have corresponding fasta file co-located)

    --index                (-idx)     Index of the job (to override LSB_JOBINDEX as used on LSF farms)

    --verbose

  OPTIONAL:

    --sampleToIgnoreInUnmatched    (-sp) Unmatched normal to ignore (to be used if the sample is one of those with a normal in the panel).

    --processid                    (-p)  Id anaylsis process to be added at a CGP specific header.

    --flagConfig                   (-c)  Config ini file to use for flag list and settings.

    --flagToVcfConfig             (-v)  Config::Inifiles style config file containing VCF flag code to flag name conversions see
                                         ../config/flag.to.vcf.convert.ini for example

    --studyType            (-t)       Study type, used to decide parameters in file (genome|genomic|WGS|pulldown|exome|WXS|followup|AMPLICON|targeted|RNA_seq).

    --flag-mnv             (-r)       Include flagging of MNVs (currently disabled)

  Examples:

    cgpFlagCaVEMan.pl [-h] -f vcfToFlag.vcf -o flaggedVCF.vcf -c configFile.ini -s human -t pulldown


=cut
