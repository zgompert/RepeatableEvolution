#!/usr/bin/perl
#
# fit gemma BSLMM for M. sativa, non-chemistry, real and randomized 
#

use Parallel::ForkManager;
my $max = 80;
my $pm = Parallel::ForkManager->new($max);

$g = "fha2013_gen_full_gus.geno";
#$g = "fha2013_gen_gus.geno";

foreach $p (@ARGV){
	foreach $ph (1..1){ 
	#foreach $ph (1..3){ 
		foreach $ch (0..9){
			sleep 2;
			$pm->start and next;
			if($ch==0){
				$o = "o_fha2013_mel_gus_lmm_ph$ph";
				#$o = "o_fha2013_gus_lmm_ph$ph";
    				system "gemma -g $g -p $p -gk 1 -o $o -maf 0\n";
    				system "gemma -g $g -p $p -k output/$o.cXX.txt -lmm 4 -n $ph -o $o -maf 0\n";
			}
			$o = "o_fha2013_mel_gus_bslmm_ph$ph"."_ch$ch";
			#$o = "o_fha2013_gus_bslmm_ph$ph"."_ch$ch";
    			system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 200000 -s 1000000\n";
			$pm->finish;
		}
	}
}
$pm->wait_all_children;

