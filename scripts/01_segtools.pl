#!/usr/bin/perl
opendir DIR,"/project/zhangq/fuxin/SEGtool_result/single_sample_result"||die"Can't open segtool result dir!";
@seg_path=readdir DIR;

foreach $dir(@seg_path){
	if($dir=~/^\./){next;}
		else{
		#	print"It's $dir thing~";
			open FILE,"</project/zhangq/fuxin/SEGtool_result/single_sample_result/$dir";
			open  OUT,">/home/fux/fux/github/ImmuneTCGA/expression_analys/tissue.specif.tcga/segtool_result_high/$dir";
			while(<FILE>){
				if(/TS_low_gene/){
		#	print"Cut $dir";
				last;}
				else{
					print OUT"$_";
					}
			}
		close FILE;
		close OUT;

#			print "$_";
#			print "$dir";
		}
}

close DIR;


