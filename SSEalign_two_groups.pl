#!/usr/bin/perl
# Description: This script is for Secondary Structure Element Alignment for two groups of proteins
# usage example: perl SSEalign_two_groups.pl <Query_SSE.fasta> <Subject_SSE.fasta> <output.txt> <Widen_cutoff(suggest 59.18)> <Length_ratio_cutoff(suggest 2) 
# Author: YANG Zhiyuan ( yangzhiyuan@link.cuhk.edu.hk )
# Institute: The Chinese University of Hong Kong
# Date: 2017-04-12

use strict;
use warnings;

# check the input parameters
 if ( $#ARGV<2 )  {   die "\n Please input 5 parameters as following:\nperl SSEalign_two_groups.pl <Query_SSE.fasta> <Subject_SSE.fasta> <output.txt> [Widen_cutoff(optional)] \n\n Other suggestions:\n1.The <Query_SSE.fasta> and <Subject_SSE.fasta> must be in FASTA format;\n2.The names of the sequences in the files <Query_SSE.fasta> and <Subject_SSE.fasta> must be only composed of English letters,numbers and underlines; \n3. The parameter [Widen_cutoff] is optional and default value is setting as 59.18 because the FDR vaule is 0.01 when Widen_cutoff is 59.18%. \n\n";    }

 unless (-d 'tmp' ) {  system ('mkdir tmp');   }
 
 # check the files
 unless (-e $ARGV[0] ) { die "The file $ARGV[0] is not existent! \n";  }
 unless (-e $ARGV[1] ) { die "The file $ARGV[1] is not existent! \n";  } 
 unless (-d 'tmp' ) {  system ('mkdir tmp');   }
 
 my $QUERY=$ARGV[0];
 my $SUBJECT=$ARGV[1];
 my $OUTFILE=$ARGV[2];
 
 my $Widen_cutoff;
 if ( $#ARGV==2 )  {  
    $Widen_cutoff=59.18;  
 }   elsif  ( $ARGV[3]<0 or $ARGV[3]>100 )  {
    die "The Widen_cutoff must be in the range [0,100]! \n";
 }   else  {
    $Widen_cutoff=$ARGV[3];
 }
  

 my %seq1;
 my %seq2;

open (OUT,">$OUTFILE") or die "cannot open $OUTFILE \n";

print OUT "Query\tSubject\tGap number\tGap open number\tPercent\tWiden%\n";

{
local $/="\n>";
open (QUE,"<$QUERY"); 
 while (my $content=<QUE>)  {
  $content=~s/>//g;
  my $loc=index($content,"\n",0);
  my $firstline=substr($content,0,$loc);
  my @parts=split/\s+/,$firstline;
  
  my $sequence=substr($content,$loc+1);
  $sequence=~s/\n//g;
  $sequence=uc($sequence);
  $seq1{$parts[0]}=$content;
} 
 close QUE;
 
 open (SUB,"<$SUBJECT") or die "cannot open $SUBJECT \n";
 while (my $content=<SUB>)  {
  $content=~s/>//g;
  my $loc=index($content,"\n",0);
  my $firstline=substr($content,0,$loc);
  my @parts=split/\s+/,$firstline;
  
  my $sequence=substr($content,$loc+1);
  $sequence=~s/\n//g;
  $sequence=uc($sequence);
  $seq2{$parts[0]}=$content;
} 
 close SUB;
}

my $suc_count=0;
my $pro_count=0;
my @querys=sort keys %seq1;
my @subjects=sort keys %seq2;

 foreach my $query(@querys) {
   foreach my $subject(@subjects) {
      $pro_count++;
      my $name="$query\_$subject";
          $name=~s/\W//g;
      my $queseq="";;
      my $hitseq="";
      my $align_seq="";
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

       open (TEMP,">tmp/$name.fasta") or die "cannot open $name.fasta \n";
          print TEMP ">$seq1{$query}\n"; 
          print TEMP ">$seq2{$subject}\n"; 
       close TEMP;

       $align_seq=`./pkg/seq-align/bin/needleman_wunsch --file tmp/$name.fasta --substitution_matrix ./pkg/seq-align/scoring/SSEalign_Matrix.txt --gapopen -2737 --gapextend -1157`;
         
       if ( $align_seq ne "" )  {

         my @seqs=split/\n/,$align_seq;
         $queseq=$seqs[0];
         $hitseq=$seqs[1]; 
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
             $suc_count++;  
             print OUT "$query\t$subject\t$gaps\t$opennum\t$percent\t$Widen\n";   
          }
        }
    }  
  }

print "\n $suc_count out of total $pro_count SSE sequence alignments pass the cutoff value!\n";  
  
 close OUT;
 
