#!/usr/local/bin/perl

#use strict;
#use warnings;

($#ARGV==4) || die"usage: [DataSetName subdirname globalTag laserTag runList(132440,132442,..)]\n";


my $dataset = $ARGV[0];
my $datasetname = $ARGV[1];
my $GRTagName = $ARGV[2];
my $laserTagName = $ARGV[3];
my $runlist = $ARGV[4];

my @lists = split(/,/, $runlist);

system("setenv SCRAM_ARCH slc5_amd64_gcc434");

if( $laserTagName == "-"){
    print "running without additional laser Tag only $GRTagName \n";
}else{
    print "running with laserTag $laserTagName and $GRTagName \n"
}

foreach (@lists) {
    my $run = $_; 
    print $run, "\n";
    
    my $dir1 = substr $run,0,3; 
    my $dir2 = substr $run,3,3; 
    print "checking ",$dir1,$dir2,"\n" ; 


    my $newdir2 = "crab_alca/${datasetname}";
    if( -d $newdir2){
    }else{
	system("mkdir $newdir2");
    }
    my $newdir = "crab_alca/${datasetname}/$run" ;
    if (-d $newdir) {
	print "directory already created.. $newdir \n";
	#system("/bin/rm -rf $newdir");
	next; 
    }
    print "making directory $newdir..\n"; 
    system("mkdir $newdir");

    if( $laserTagName == "-"){
	system("cp crab_aclcapiz_temp/runpizfromdatarawnolasertag.py $newdir/runpizfromdata.py");
	system("cp crab_aclcapiz_temp/crab.cfg.lpcraw2reco $newdir/crab.cfg");
    }else{
	system("cp crab_aclcapiz_temp/crab.cfg.lpcraw2recolasertag $newdir/crab.cfg");
	system("cp crab_aclcapiz_temp/runpizfromdataraw.py $newdir/runpizfromdata.py");
    }

    my $outdir = "/pnfs/cms/WAX/resilient/username/data/pizdata/$run";
    if(-d $outdir) {
        print "directory in resilient already exist..\n";
        system("chmod 775 $outdir");
    }else{
        print "making directory in resilient ..\n";
        system("mkdir $outdir");
        system("chmod 775 $outdir");
    }
    
    system("perl ReplaceString.pl InputDataSet $dataset $newdir/crab.cfg");
    system("perl ReplaceString.pl datasetname $datasetname $newdir/crab.cfg");
    system("perl ReplaceString.pl xxxxxx $run $newdir/crab.cfg");
    system("perl ReplaceString.pl xxxxxx $run $newdir/runpizfromdata.py");
    system("perl ReplaceString.pl datasetname $datasetname $newdir/runpizfromdata.py");
    system("perl ReplaceString.pl laserTagName $laserTagName $newdir/runpizfromdata.py");
    system("perl ReplaceString.pl GRTagName $GRTagName $newdir/runpizfromdata.py");
    system("perl ReplaceString.pl laserTagName $laserTagName $newdir/crab.cfg");
    system("perl ReplaceString.pl GRTagName $GRTagName $newdir/crab.cfg");
    
    print "crab directory ready...\n"; 
    
    system("cd $newdir; pwd; source /uscmst1/prod/grid/gLite_SL5_CRAB_27x.sh; eval `scramv1 runtime -sh`; which cmsRun;source /uscmst1/prod/grid/CRAB/crab.sh; which crab; crab -create -submit");
    
    print "crab job submitted for run $run ..\n"; 
    
    
}
