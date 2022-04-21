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

package Sanger::CGP::CavemanPostProcessor::ConfigParser;

use Sanger::CGP::CavemanPostProcessor::Constants;

use YAML::XS qw(LoadFile DumpFile);
use Config::IniFiles;
use Const::Fast qw(const);
use File::Basename qw(fileparse);
use List::Util qw(first);
use Data::Dumper;
use Carp;

my $const = 'Sanger::CGP::CavemanPostProcessor::Constants';

sub new {
  my ($proto) = shift;
  my %inputs = @_;
  my $class = ref($proto) || $proto;
   my $self = {};
  bless($self, $class);
  return $self;
}

sub _get_cfg_ini_config_params{
  my ($file, $species, $seq_type, $flagopts) = @_;
  my $cfg = Config::IniFiles->new( -file => $file, -allowcontinue => 1);
  #Get parameter group.
  my $alternate = "";
  my $sppTypeCombo = "".$species."_".$seq_type;
  my $paramSectName = $sppTypeCombo." ".$const->config_param_names('CONFIG_PARAMETERS');
  my @parameterNames = $cfg->Parameters($paramSectName);
  #iterate through each parameter name and load into a hash for module setup.
  my ($sectParams,$bedFileParams);
  foreach my $paramName(@parameterNames){
    $sectParams->{$paramName} = $cfg->val($paramSectName,$paramName);
  }
  #Get the flaglist group
  my @flagList;
  if(! defined($flagopts)){
    #Get the flaglist group
    $paramSectName = $sppTypeCombo." ".$const->config_param_names('CONFIG_FLAGLIST');
    my $prm = $const->config_param_names('FLAG_PARAMETER');
    @flagList = $cfg->val($paramSectName,$prm);
    if(!@flagList){
      croak("No flagList found in ".$file." for section $paramSectName. No flagging will be done.");
      @flagList = ();
      return ($sectParams,\@flagList,$bedFileParams,undef);
    }
  }else{
    @flagList = @$flagopts;
  }
  if(!defined($sectParams)){
    croak("No config found in ".$opts->{'c'}." for section $paramSectName");
  }
  #Get the bedfiles
  $paramSectName = $sppTypeCombo." ".$const->config_param_names('CONFIG_BEDFILES');
  @parameterNames = $cfg->Parameters($paramSectName);
  foreach my $pName(@parameterNames){
    $bedFileParams->{$pName} = $cfg->val($paramSectName,$pName);
  }
  if(!defined($bedFileParams)){
    croak("No bed file parameters found in ".$opts->{'c'}." for section $paramSectName");
  }
  #Get the MNV flaglist group. This COULD be non existant
  my @mnv_flaglist;
  $paramSectName = $sppTypeCombo." ".$const->config_param_names('CONFIG_MNVFLAGLIST');
  $prm = $const->config_param_names('FLAG_PARAMETER');
  @mnv_flaglist = $cfg->val($paramSectName,$prm);
  if(!@mnv_flaglist){
    #carp("No ".$const->config_param_names('CONFIG_MNVFLAGLIST')." found in ".$file." for section $paramSectName. No flagging will be done.");
    @mnv_flaglist = ();
    return ($sectParams,\@flagList,$bedFileParams,\@mnv_flaglist);
  }
  return ($sectParams, \@flagList, $bedFileParams, \@mnv_flaglist);
}

sub _get_yml_config_params{
  my ($file, $species, $seq_type) = @_;
  my $hashref_yaml = LoadFile( $file );
  my $section = $hashref_yaml->{$species}->{$seq_type};
  if(! $section){
    croak("No section combination for $species $seq_type found in config file '$file'.");
  }
  my $sectParams = $section->{$const->config_param_names('CONFIG_PARAMETERS')};
  if(! $sectParams){
    croak("No ".$const->config_param_names('CONFIG_PARAMETERS')." section for $species $seq_type found in config file '$file'.");
  }
  my $flagList = $section->{$const->config_param_names('CONFIG_FLAGLIST')};
  if(! $flagList){
    croak("No ".$const->config_param_names('CONFIG_FLAGLIST')." section for $species $seq_type found in config file '$file'.");
  }
  # MNV may not exist so only warn
  my $mnvflaglist = $section->{$const->config_param_names('CONFIG_MNVFLAGLIST')};
  if(! $mnvflaglist){
    carp("No ".$const->config_param_names('CONFIG_MNVFLAGLIST')." section for $species $seq_type found in config file '$file'.");
  }
  my $bedFileParams = $section->{$const->config_param_names('CONFIG_BEDFILES')};
  if(! $bedFileParams){
    croak("No ".$const->config_param_names('CONFIG_BEDFILES')." section for $species $seq_type found in config file '$file'.");
  }
  return ($sectParams, $flagList, $bedFileParams, $mnvflaglist);
}

sub getConfigParams{
  my ($opts, $flagopts) = @_;
  #Determine input config type
  my($base, $dir, $ext) = fileparse($opts->{'c'}, qw{ .ini .yaml .yml});
  if($ext eq '.yaml' || $ext eq '.yml'){ # Just in case we encounter a .yml rather than yaml suffix
    return _get_yml_config_params($opts->{'c'},$opts->{'s'},$opts->{'t'});
  }elsif($ext eq '.ini'){
    return _get_cfg_ini_config_params($opts->{'c'},$opts->{'s'},$opts->{'t'}, $flagopts);
  }else{
    croak("Incorrect file extension for ini file '".$opts->{'c'}."'. Expected .yaml or .ini but got '$ext'");
  }
}

sub validate_yaml_config{
  my ($input_file) = @_;
  # Validate yaml input file format
  my $hashref_yaml = LoadFile( $file );
  my $section = $hashref_yaml->{$species}->{$seq_type};
  if(! $section){
    croak("No section combination for $species $seq_type found in config file '$file'.");
  }
  my $sectParams = $section->{$const->config_param_names('CONFIG_PARAMETERS')};
  if(! $sectParams){
    croak("No ".$const->config_param_names('CONFIG_PARAMETERS')." section for $species $seq_type found in config file '$file'.");
  }
  my $flagList = $section->{$const->config_param_names('CONFIG_FLAGLIST')};
  if(! $flagList){
    croak("No ".$const->config_param_names('CONFIG_FLAGLIST')." section for $species $seq_type found in config file '$file'.");
  }
  my $bedFileParams = $section->{$const->config_param_names('CONFIG_BEDFILES')};
  if(! $bedFileParams){
    croak("No ".$const->config_param_names('CONFIG_BEDFILES')." section for $species $seq_type found in config file '$file'.");
  }
  return
}

sub convert_ini_to_yaml{
  #Convert old ini file to yaml format
  my ($input_ini, $output_yaml_file) = @_;
  #Read in ini file into hash
  my $cfg = Config::IniFiles->new( -file => $input_ini, -allowcontinue => 1);
  my %species_seq_type;
  foreach my $section ($cfg->Sections){
    my ($species_type, undef) = split /\s+/ , $section ; 
    my ($spp, $seq_type) = split /_/ , $species_type;
    if (! exists($species_seq_type{$spp})){
      $species_seq_type{$spp} = []
    }
    push(@{$species_seq_type{$spp}},$seq_type) unless (defined( first {$_ eq $seq_type} @{$species_seq_type{$spp}}));
  }
  my $yaml_out;
  foreach my $species(keys(%species_seq_type)){
    foreach my $seq_type(@{$species_seq_type{$species}}){
      my $key = $species."_".$seq_type;
      $yaml_out->{$species}->{$seq_type} = {};
      my $paramSectName = $key." ".$const->config_param_names('CONFIG_PARAMETERS');
      my @parameterNames = $cfg->Parameters($paramSectName);
      #iterate through each parameter name and load into a hash for module setup.
      my ($sectParams,$bedFileParams);
      foreach my $paramName(@parameterNames){
        $yaml_out->{$species}->{$seq_type}->{$const->config_param_names('CONFIG_PARAMETERS')}->{$paramName} = $cfg->val($paramSectName,$paramName);
      }

      $paramSectName = $key." ".$const->config_param_names('CONFIG_FLAGLIST');
      my $flagprm = $const->config_param_names('FLAG_PARAMETER');
      my @flagList = $cfg->val($paramSectName,$flagprm);
      $yaml_out->{$species}->{$seq_type}->{$const->config_param_names('CONFIG_FLAGLIST')} = \@flagList;

      $paramSectName = $key." ".$const->config_param_names('CONFIG_MNVFLAGLIST');
      my $flagmnvprm = $const->config_param_names('FLAG_PARAMETER');
      if ($cfg->exists($paramSectName, $flagmnvprm)){
        my @mnvflagList = $cfg->val($paramSectName,$flagmnvprm);
        $yaml_out->{$species}->{$seq_type}->{$const->config_param_names('CONFIG_MNVFLAGLIST')} = \@mnvflagList;
      }

      $paramSectName = $key." ".$const->config_param_names('CONFIG_BEDFILES');
      @parameterNames = $cfg->Parameters($paramSectName);
      foreach my $pName(@parameterNames){
        $val = $cfg->val($paramSectName,$pName);
        if ($val eq ''){
          $val = undef 
        }
        $yaml_out->{$species}->{$seq_type}->{$const->config_param_names('CONFIG_BEDFILES')}->{$pName} = $val;
      }
    }
  }
  #Print yaml to file
  DumpFile($output_yaml_file,$yaml_out);
  return;
}

return 1;
