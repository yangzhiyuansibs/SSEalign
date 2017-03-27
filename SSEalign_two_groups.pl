#!/usr/bin/perl
# Description: This script is for Secondary Structure Element Alignment for two groups of proteins
# usage example: perl SSEalign_two_groups.pl <Query_SSE.fasta> <Subject_SSE.fasta> <output.txt> <Widen_cutoff(suggest 59.18)> <Length_ratio_cutoff(suggest 2) 
# Author: YANG Zhiyuan ( yangzhiyuan@link.cuhk.edu.hk )
# Institute: The Chinese University of Hong Kong
# Date: 2017-03-23

use strict;
use warnings;

# check the input parameter
 if ( $#ARGV<4 )  {   die "\n Please input 5 parameter as following:\nperl SSEalign_two_groups.pl <Query_SSE.fasta> <Subject_SSE.fasta> <output.txt> <Widen_cutoff(suggest 59.18)> <Length_ratio_cutoff(suggest 2)> \n\n Other suggestions:\n1.The <Query_SSE.fasta> and <Subject_SSE.fasta> must be in FASTA format;\n2.The name contained NON-word character such as sp|P0ABI8|CYOB_ECOLI is not allowed in <Query_SSE.fasta> and <Subject_SSE.fasta>; \n3.The tmp folder is also needed; \n4.The FDR is 0.01 when Widen_cutoff is 59.18%, so this value was suggested. \n5 Length_ratio_cutoff was calculated by the longer length divide the shorter length,so it must larger than 1. \n\n";    }

 if ( $ARGV[3]<1 or $ARGV[3]>100 )  { die "The Widen_cutoff must be in the range [1,100]! \n"; }
 if ( $ARGV[4]<1 or $ARGV[4]>20  )  { die "The Length_ratio_cutoff must be in the range [1,20]! \n"; }
 unless (-d 'tmp' ) {  system ('mkdir tmp');   }
 
 
 # check the file
 unless (-e $ARGV[0] ) { die "The file $ARGV[0] is not existent! \n";  }
 unless (-e $ARGV[1] ) { die "The file $ARGV[1] is not existent! \n";  } 
 unless (-d 'tmp' ) {  system ('mkdir tmp');   }
 
 my $QUERY=$ARGV[0];
 my $SUBJECT=$ARGV[1];
 my $OUTFILE=$ARGV[2];
 my $Widen_cutoff=$ARGV[3];
 my $Len_cutoff=$ARGV[4];
 
 my @querys;
 my @subjects;
 my %len;
 my %sequ;

open (OUT,">$OUTFILE") or die "cannot open $OUTFILE \n";

print OUT "Query\tSubject\tHH\tEE\tCC\tHE\tHC\tEC\tGaps\tGap_open\tIdentity\tScore\n";

{
local $/="\n>";
open (QUE,"<$QUERY"); 
 while (my $content=<QUE>)  {
  $content=~s/>//g;
  my $loc=index($content,"\n",0);
  my $firstline=substr($content,0,$loc-1);
  my @parts=split/\s+/,$firstline;
  
  my $sequence=substr($content,$loc+1);
  $sequence=~s/\n//g;
  $sequence=uc($sequence);
  $len{$parts[0]}=length($sequence);
  open (TEMP,">tmp/$parts[0].fasta");
      print TEMP ">Query.$content"; 
  close TEMP;
} 
 close QUE;
 
 open (SUB,"<$SUBJECT") or die "cannot open $SUBJECT \n";
 while (my $content=<SUB>)  {
  $content=~s/>//g;
  my $loc=index($content,"\n",0);
  my $firstline=substr($content,0,$loc-1);
  my @parts=split/\s+/,$firstline;
  
  my $sequence=substr($content,$loc+1);
  $sequence=~s/\n//g;
  $sequence=uc($sequence);
  $len{$parts[0]}=length($sequence);
  open (TEMP,">tmp/$parts[0].fasta");
      print TEMP ">Subje.$content"; 
  close TEMP;

} 
 close SUB;
}

 foreach my $query(@querys) {
   print "$query\n";
   foreach my $subject(@subjects) {
      my $name="$query.$subject";
      my $queseq='';;
      my $hitseq='';
      my $percent=0;
      my $Widen=0;
      my $CC=0;
      my $HC=0;
      my $EC=0;
      my $HH=0;
      my $HE=0;
      my $EE=0;
      my $gaps=0;
      my $opennum=0;
      my @quechars;
      my @hitchars;
      unless (  $len{$query}/$len{$subject}>$Len_cutoff  or  $len{$subject}/$len{$query}>$Len_cutoff )  {

    `stretcher -asequence tmp/$query.fasta -bsequence tmp/$subject.fasta -outfile tmp/$name.stre -auto Y;`;
         
    if (-e "tmp/$name.stre" )  {
      open (CAL,"<tmp/$name.stre"); 
       while (my $line=<CAL>) {
        if ( $line=~/^Query\.\s+(\S+)\n/ ) {  $queseq=$queseq.$1;   }
        if ( $line=~/^Subje\.\s+(\S+)\n/ ) {  $hitseq=$hitseq.$1;   } 
       }
      close CAL;
  
      my @quegap=split/\-+/,$queseq;
      my @hitgap=split/\-+/,$hitseq;
      $opennum=$#quegap+$#hitgap;
      if ( substr($queseq,0,1) eq "-" )  { $opennum=$opennum-1;   }
      if ( substr($hitseq,0,1) eq "-" )  { $opennum=$opennum-1;   }   
       
      @quechars=split(//, $queseq);
      @hitchars=split(//, $hitseq);  
      my $i;
      for ($i=0; $i<$#quechars; $i++ ) {
 
        if (($quechars[$i] eq "-") or ($hitchars[$i] eq "-") ) {  $gaps++;   }

        if (($quechars[$i] eq "C") and ($hitchars[$i] eq "H") ) {  $HC++;   }
        if (($quechars[$i] eq "C") and ($hitchars[$i] eq "E") ) {  $EC++;   }
        if (($quechars[$i] eq "C") and ($hitchars[$i] eq "C") ) {  $CC++;   }

        if (($quechars[$i] eq "H") and ($hitchars[$i] eq "H") ) {  $HH++;   }
        if (($quechars[$i] eq "H") and ($hitchars[$i] eq "E") ) {  $HE++;   }
        if (($quechars[$i] eq "H") and ($hitchars[$i] eq "C") ) {  $HC++;   }

        if (($quechars[$i] eq "E") and ($hitchars[$i] eq "H") ) {  $HE++;   }
        if (($quechars[$i] eq "E") and ($hitchars[$i] eq "E") ) {  $EE++;   }
        if (($quechars[$i] eq "E") and ($hitchars[$i] eq "C") ) {  $EC++;   }   
      }

     $percent=int(($HH+$EE+$CC)/($HH+$EE+$CC+$HE+$HC+$EC+$gaps)*10000+0.5)/10000;
     $Widen=int((1.796*$HH+2.374*$EE+0.351*$CC)/(1.796*$HH+2.374*$EE+0.351*$CC+2.478*$HE+0.557*$HC+0.460*$EC+1.157*($gaps-$opennum)+2.737*$opennum)*10000+0.5)/100;

        if ( $Widen>=$Widen_cutoff )   { 
           print OUT "$query\t$subject\t$HH\t$EE\t$CC\t$HE\t$HC\t$EC\t$gaps\t$opennum\t$percent\t$Widen\n";   

           my $description="Query=$query  Subject=$subject  Gapnum=$gaps  Gapopen=$opennum  percent=$percent  Widen=$Widen";  
           `sed -i \"32i $description\" tmp/$name.stre`;
          }
          else {  `rm tmp/$name.stre`;    }
      }
    }  
  }
}

  close OUT;