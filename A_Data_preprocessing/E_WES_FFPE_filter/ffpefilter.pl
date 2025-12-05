
## pipeline for FFPE filtering ##

#!/usr/bin/perl

use strict;
use warnings;
#use POSIX;
use Getopt::Long;

my $version = 0.0.1;

my $green = "\e[32m";
my $yellow = "\e[33m";
my $normal = "\e[0m";
my $red = "\e[31m";

(my $usage = <<OUT) =~ s/\t+//g;
FFPE filtering pipeline 
Pipeline version: $version
$yellow     Usage: perl $0  --rdir --groupname --users --step 

<rdir> = full path of the folder holding files for this run (user must provide)
<log> = full path of the folder for saving log file; usually upper folder of rdir
<groupname> = job group name
<users> = user name for job group
<step> run this pipeline step by step. (user must provide)

$green       [1]  Run MicroSEC

OUT

my $step_number = -1;
my $compute_username="";
my $group_name="";
my $help = 0;
my $q_name="";
my $run_dir="";
my $log_dir="";

my $status = &GetOptions (
      "step=i" => \$step_number,
      "groupname=s" => \$group_name,
      "users=s" => \$compute_username,	
      "rdir=s" => \$run_dir,
      "q=s" => \$q_name,
      "log=s"  => \$log_dir,
      "help" => \$help, 
	);

if ($help || $run_dir eq "" || $log_dir eq "" || $group_name eq "" || $compute_username eq "" || $step_number<0) {
	 print "wrong option\n";
	  print $usage;
      exit;
   }

print "run dir=",$run_dir,"\n";
print "log dir=",$log_dir,"\n";
print "step num=",$step_number,"\n";
print "job group=",$group_name,"\n";
print "user group=",$compute_username,"\n";
print "queue name=",$q_name,"\n";

if($q_name eq "") 
{
	$q_name="dinglab";
}

my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-1];
my $HOME1=$log_dir;

if (! -d $HOME1)
{
`mkdir $HOME1`; 
}
if (! -d $HOME1."/tmpffpe") {
    `mkdir $HOME1"/tmpffpe"`;
}
my $job_files_dir = $HOME1."/tmpffpe";
#store SGE output and error files here

if (! -d $HOME1."/LSF_DIR_FFPE") {
    `mkdir $HOME1"/LSF_DIR_FFPE"`;
}
my $lsf_file_dir = $HOME1."/LSF_DIR_FFPE";


my $run_script_path =`echo \$PWD`;
chomp $run_script_path;
my $script_dir=$run_script_path; 
print $script_dir,"\n";
$run_script_path = "/usr/bin/perl ".$run_script_path."/";


print $run_script_path,"\n";
my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

### run ffpe for each sample in sample_dir_list directory

if ($step_number < 2 && $step_number>0) {
    for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
	    if($step_number == 1) {
                    &bsub_microsec(1);	
			}
			}
	}
	}
}

sub bsub_microsec{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

   $current_job_file = "j1_microsec".$sample_name.".sh";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

   my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";


    if (! -s $IN_bam_T) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_T is empty!", $normal, "\n\n";
    }

    open(MICROSEC, ">$job_files_dir/$current_job_file") or die $!;
    print MICROSEC "#!/bin/bash\n";
    print MICROSEC "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print MICROSEC "MAF=".$sample_full_path."/".$sample_name.".formated.maf\n";
    print MICROSEC "MUTCONVERTED=".$sample_full_path."/".$sample_name."_Mutations_converted.xlsx\n";
    print MICROSEC "RUNDIR=".$sample_full_path."\n";
    print MICROSEC "Rscript ".$script_dir."/convert_docker.R \$RUNDIR \$MAF $sample_name","\n";
    print MICROSEC "Rscript ".$script_dir."/MicroSEC_docker.R $sample_name $sample_name \$MUTCONVERTED \$TBAM NULL 151 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT hg38 \$RUNDIR","\n";  
    close MICROSEC;

    my $sh_file=$job_files_dir."/".$current_job_file;  
    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/microsec:0.0.1)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);
 
}

