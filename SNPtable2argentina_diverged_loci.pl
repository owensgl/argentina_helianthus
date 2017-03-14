#!/bin/perl
use warnings;
use strict;

my $popfile = $ARGV[0];

my %pop;
open POP, $popfile;
while(<POP>){
    chomp;
    my @a = split(/\t/,$_);
    $pop{$a[0]} = $a[1];
}
close POP;

my @samples;
my %names;
my %count;
my $sample_n;
while(<STDIN>){
    chomp;
    my @a = split(/\t/,$_);
    my %a_freq;
    my %calls;
    my %total_alleles;
    if ($. == 1){
        foreach my $i (2..$#a){
            push (@samples, $a[$i]);
            $names{$i} = $a[$i];
	    $sample_n = $#a;
        }
    }else{

        foreach my $i (2..$#a){
            if (($pop{$names{$i}}) and ($a[$i] ne "NN")){
                my @bases = split(//, $a[$i]);
                $a_freq{$pop{$names{$i}}}{$bases[0]}++;
                $a_freq{$pop{$names{$i}}}{$bases[1]}++;
                $calls{$pop{$names{$i}}}+=2;
		$total_alleles{$bases[0]}++;
		$total_alleles{$bases[1]}++;
            }
        }
	if (keys %total_alleles ne 2){
		next;
	}
        my @P1_bases = sort { $a_freq{"P1"}{$b} <=> $a_freq{"P1"}{$a} } keys %{$a_freq{"P1"}} ;
        my $P1_major = $P1_bases[0];
#	print "$P1_bases[0]\t$P1_bases[1]\n";
        my @P2_bases = sort { $a_freq{"P2"}{$b} <=> $a_freq{"P2"}{$a} } keys %{$a_freq{"P2"}} ;
        my $P2_major = $P2_bases[0];
        my $P1_freq = $a_freq{"P1"}{$P1_major}/ $calls{"P1"};
        my $P2_freq = 1- $a_freq{"P2"}{$P2_major}/ $calls{"P2"};
        #print "$P1_freq:$P2_freq\n";
	foreach my $i (2..$#a){
	    if ($a[$i] ne "NN"){
		my @bases = split(//, $a[$i]);
		foreach my $n (0..1){
		    if ($bases[$n] eq $P1_major){
			$count{$names{$i}}{"P1"}++;
		    }elsif($bases[$n] eq $P2_major){
			$count{$names{$i}}{"P2"}++;
		    }
		}
	    }
	}
    }
}

print "name\tP1\tP2\ttotal\tpercent";
foreach my $i (2..$sample_n){
	unless ($count{$names{$i}}{"P1"}){
        	$count{$names{$i}}{"P1"} = 0;
        }
        unless ($count{$names{$i}}{"P2"}){
                $count{$names{$i}}{"P2"} = 0;
        }
        my $total = $count{$names{$i}}{"P1"} + $count{$names{$i}}{"P2"};
	my $percent = $count{$names{$i}}{"P1"}/ $total;
        print "\n$names{$i}\t$count{$names{$i}}{'P1'}\t$count{$names{$i}}{'P2'}\t$total\t$percent";
}
