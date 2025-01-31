#!/usr/bin/perl

$fn1 = $fn2 = $out1 = $out2 = $ARGV[1];
$fn1="$fn1".".passed.fastq";
$fn2="$fn2".".failed.fastq";
$out1="$out1".".passed.out";
$out2="$out2".".failed.out";
#$fn1 =~ s/(\S+).fastq$/\1.passed.fastq/g;
#$fn2 =~ s/(\S+).fastq$/\1.failed.fastq/g;
#$out1 =~ s/(\S+).fastq$/\1.passed.out/g;
#$out2 =~ s/(\S+).fastq$/\1.failed.out/g;

print "$ARGV[0]\t$ARGV[1]\t$fn1\n";

open (my $passed, ">", $fn1 ) or die "Could not open file $fn1";
open (my $failed, ">", $fn2 ) or die "Could not open file $fn2";
open (my $out1, ">", $out1 ) or die "Could not open file $out1";
open (my $out2, ">", $out2 ) or die "Could not open file $out2";


while (<>) {
  chomp();
  if (/^(@\S+)\s\S+$/){
    $rid=$1;
    $inSeq=1;
    next;
  }
#  if ($inSeq!=1) { next;}
  if ($inSeq==1) {
  $pcr1=$link1=$seq1=$rem1="NA";
  $M1=0;
  ## match guide RNA 1.
  if (/^([ACTGN]{0,15})(TTGTGGAAAGGACGAAACACCG)([ACTGN]{18,21})(GT\S+)/) {

    ($pcr1,$link1,$seq1,$rem1,$pos1,$pos2,$pos3,$pos4)=($1,$2,$3,$4,$-[1],$+[3]);
    $M1 = 1;
  } else { $S1 = $_; }
  }
  if ($inSeq==3){
    if ( $M1 == 1 ){
      $quality = substr $_, $pos1, 19;
      #$quality2= substr $_, $pos3, $pos4-$pos3;
      print $out1 "$rid\t$pcr1\t$link1\t$seq1\t$rem1\t$pos1\t$pos2\n";#\t$quality\n";
      #print $passed "$rid\n${pcr1}$link1$seq1\n+\n$quality\n"
      print $passed "$rid\n$seq1\n+\n$quality\n"
  }
    else { 
      #  print $passed "$rid BC:Z:NNNNN\nNNNNN\n+\nNNNNN\n";
      print $out2 "$rid\t$S1\t$_\n"; 
      print $failed "$rid\n$S1\n+\n$_\n";
      }
  }
$inSeq++;
}

