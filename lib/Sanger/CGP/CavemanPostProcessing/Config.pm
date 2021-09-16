##########LICENCE##########
# Copyright (c) 2014-2018 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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
package Sanger::CGP::CavemanPostProcessing::Config;

use strict;
use warnings FATAL=>'all';
use autodie;
use Const::Fast;
use Carp qw(croak cluck);
use File::Spec;
use Data::Dumper;
use Config::IniFiles;

use Sanger::CGP::CavemanPostProcessing;
use Sanger::CGP::CavemanPostProcessing::FlagConfig;

const my $CONFIG_DEFAULT => 'DEFAULT';
const my $CONFIG_FLAGLIST => "%s FLAGLIST";
const my $CONFIG_PARAMETERS => "%s PARAMS";
const my $CONFIG_BEDFILES => "%s BEDFILES";

const my $CONFIG_INFO => "info";
const my $CONFIG_ID => "id";
const my $CONFIG_DESCRIPTION => "description";
const my $CONFIG_VAL => "val";
const my $CONFIG_TYPE => "type";
const my $CONFIG_INTERSECT => "intersect";
const my $CONFIG_NEED_VCF => "needs_vcf";
const my $CONFIG_OPTIONNAME => "optname";
const my $CONFIG_FILENAME => "filename";
const my $CONFIG_NUMBER => 'Number';

const my $VCF_UM_FLAG_NAME => 'unmatchedNormalVcfFlag';

sub new {
	my ($proto) = @_;
	my $class = ref($proto) || $proto;
 	my $self = {};
  bless($self, $class);
  return $self;
}

sub init_config{
  my ($self,$opts) = @_;
  return $self->_setup($opts);
}

sub _setup{
  my ($self,$opts) = @_;
  my $cfg = Config::IniFiles->new( -file => $opts->{'c'}, -allowcontinue => 1, -default => $CONFIG_DEFAULT );
  $self->{'cfg'} = $cfg;
  #Check config contains a group named appropriately
  croak("Species or type not defined in optionts hash.") unless(exists($opts->{'s'}) && exists($opts->{'t'}) && defined($opts->{'s'}) && defined($opts->{'t'}));
  my $group = $opts->{'s'}."_".$opts->{'t'};
  croak("Group $group not found in config file $opts->{'c'}.") unless(grep {/$group/} $self->{'cfg'}->Groups);
  #Check for relevant sections
  my $flagList_sect = sprintf($CONFIG_FLAGLIST,$group);
  croak("Config section $flagList_sect not found in $opts->{'c'}.") unless($self->{'cfg'}->SectionExists($flagList_sect));
  $self->{'flist_sect'} = $flagList_sect;
  @{$self->{'flist'}} = $self->{'cfg'}->Parameters($flagList_sect);
  my $param_sect = sprintf($CONFIG_PARAMETERS,$group);
  croak("Config section $param_sect not found in $opts->{'c'}.") unless($self->{'cfg'}->SectionExists($param_sect));
  $self->{'params_sect'} = $param_sect;
  @{$self->{'params'}} = $self->{'cfg'}->Parameters($param_sect);
  #Get flag information from flag_to_vcf config
  my $flagcfg = Config::IniFiles->new(-file => $opts->{'v'}, -allowcontinue => 1, -default => $CONFIG_DEFAULT );
  $self->{'flagcfg'} = $flagcfg;
  $self->{'flags'} = $self->_read_flag_config($opts);
  return;
}

sub _read_flag_config{
  my ($self,$opts) = @_;
  #Get list of flags from config
  my @flag_list;
  foreach my $flag_section(@{$self->flag_list_val($self->flaglist->[0])}){
    my $flag = Sanger::CGP::CavemanPostProcessing::FlagConfig->new();
    $flag->name($flag_section);
    my $info=0;
    if($self->flag_config->exists($flag_section,$CONFIG_INFO)){
      $info = $self->flag_config->val($flag_section,$CONFIG_INFO);
    }
    $flag->is_info($info);
    if(!($self->flag_config->exists($flag_section,$CONFIG_ID))){
      croak("No flag id found for $flag_section");
    }
    my $id = $self->flag_config->val($flag_section,$CONFIG_ID);
    $flag->id($id);
    my $val = undef;
    if($self->flag_config->exists($flag_section,$CONFIG_VAL)){
      $val = $self->flag_config->val($flag_section,$CONFIG_VAL);
    }
    $flag->value($val);
    my $type = undef;
    if($self->flag_config->exists($flag_section,$CONFIG_TYPE)){
      $type = $self->flag_config->val($flag_section,$CONFIG_TYPE);
    }
    $flag->type($type);
    my $intersect = 0;
    if($self->flag_config->exists($flag_section,$CONFIG_INTERSECT)){
      $intersect = $self->flag_config->val($flag_section,$CONFIG_INTERSECT);
      #This should then specify the option name and (optional) file path to find these files.
      croak("Option name parameter not provided for intersect flag '$flag_section'.") unless($self->flag_config->exists($flag_section,$CONFIG_OPTIONNAME));
      my $opt = $self->flag_config->val($flag_section,$CONFIG_OPTIONNAME);
      my $fname = undef;
      $fname = $self->flag_config->val($flag_section,$CONFIG_FILENAME) if($self->flag_config->exists($flag_section,$CONFIG_FILENAME));
      my $file = $opts->{$opt};
      if(defined ($fname)){
        $file = File::Spec->catfile($opts->{$opt},$fname);
      }
      #Check this has the relevant file available if it's an intersect flag.
      croak("File '$file' required by flag '$flag_section' doesn't exist or isn't readable.") unless (-e $file && -r $file);
      $flag->intersect_file($file);
    }
    $flag->is_intersect($intersect);
    if(!($self->flag_config->exists($flag_section,$CONFIG_DESCRIPTION))){
      croak("No flag description found for $flag_section");
    }
    my $need_vcf = 0;
    if($self->flag_config->exists($flag_section,$CONFIG_NEED_VCF)){
        $need_vcf = $self->flag_config->val($flag_section,$CONFIG_NEED_VCF);
    }
    $flag->is_need_vcf($need_vcf);
    my $pre_description = $self->flag_config->val($flag_section,$CONFIG_DESCRIPTION);
    #Populate the description with the config values.
    my $desc = $self->description_with_params($flag_section,$pre_description);
    $flag->description($desc);
    if($flag_section eq $VCF_UM_FLAG_NAME){
      #Special case for VCF unmatched normal, only if there's a bed file supplied. Otherwise it's old style with one VCF per contig
      if($self->flag_config->exists($flag_section,$CONFIG_FILENAME)){
        croak("Option name parameter not provided for intersect flag '$flag_section'.") unless($self->flag_config->exists($flag_section,$CONFIG_OPTIONNAME));
        my $opt = $self->flag_config->val($flag_section,$CONFIG_OPTIONNAME);
        my $fname = undef;
        $fname = $self->flag_config->val($flag_section,$CONFIG_FILENAME) if($self->flag_config->exists($flag_section,$CONFIG_FILENAME));
        my $file = $opts->{$opt};
        if(defined ($fname)){
          $file = File::Spec->catfile($opts->{$opt},$fname);
        }
        $flag->intersect_file($file);
      }
    }


    push(@flag_list,$flag);
  }#End of iteration through each requested flag
  if(!@flag_list){
		cluck("No flag list found in ".$opts->{'c'}." for section. No flagging will be done.");
		@flag_list = ();
	}
  return \@flag_list;
}

sub description_with_params{
	my ($self,$flag_section,$desc) = @_;
	foreach my $pName(@{$self->params()}){
	  if($desc =~ m/$pName/){
	    my $val = $self->param_val($pName);
      $desc =~ s/$pName/$val/g;
	  }
	}
	return $desc;
}

sub flag_config{
  my ($self) = @_;
  return $self->{'flagcfg'};
}

sub config{
  my ($self) = @_;
  return $self->{'cfg'};
}

sub flags{
  my ($self) = @_;
  return $self->{'flags'};
}

sub flaglist{
  my ($self) = @_;
  return $self->{'flist'};
}

sub params{
  my ($self) = @_;
  return $self->{'params'};
}

sub flag_list_val{
  my ($self,$param) = @_;
  my @flagList = $self->{'cfg'}->val($self->{'flist_sect'},$param);
  return \@flagList;
}

sub param_val{
  my ($self,$param) = @_;
  return $self->{'cfg'}->val($self->{'params_sect'},$param);
}

1;
