#!/usr/local/bin/perl
# perl script ReplaceString
# author: Sean Kelly UCLA - mod Athen Trattner UCB
# description: replaces old string with new string in file(s)
($#ARGV>=2) || die"usage: ReplaceString  [old string] [new string] [input
files]\n";
$oldstring = $ARGV[0];
$newstring = $ARGV[1];
for ($i=2;$i<@ARGV;$i++){
    $found=0;
    system("/bin/cp", $ARGV[$i], "$ARGV[$i].temp01");
    open(FILE,"$ARGV[$i].temp01") ||
        die "\ncan't open file $ARGV[$i].temp01\n";
    while($line = <FILE>){
        if ($line =~ (/$oldstring/)) {
            $found=1;
            print "\nfound in file $ARGV[$i]\n";
            last;
        }
    }
    close(FILE);

    if ($found) {
        open(FILE,"$ARGV[$i]") ||
            die "\ncan't open file $ARGV[$i]\n";
        open(FILE2,">$ARGV[$i].temp01.txt") ||
            die "\ncan't open file $ARGV[$i].temp01.txt\n";
        while($line = <FILE>){
            $line =~ s/$oldstring/$newstring/g;
            print FILE2 $line;
        }
        if(-e "$ARGV[$i].temp01.txt"){
            system("/bin/mv $ARGV[$i].temp01.txt $ARGV[$i]");
	}
        close(FILE2);
    }
    system("/bin/rm $ARGV[$i].temp01");
}





close(FILE);


