# Author: Hengyi Jiang <hengyi.jiang@gmail.com> <jhy@fudan.edu.cn>
# purpose: to make an index file for each fasta result, usually big fasta file, with lots of genomes together
# usage: ruby fasta_index.rb fasta_file > fasta.index
# fasta.index : each line composed of >seq_name ||start_pointer||end_pointer

if ARGV.size < 1
	puts "ruby fasta_index.rb fasta_file > fasta.index"
end

filename = ARGV.shift

begin
	f= File.open(filename);
rescue
    puts 'file open failed';
    abort
end 

started=0
name=""
pointer_left=-1;
pointer_right=-1;
seq_left=-1
width=-1
n=0
while l=f.gets do
	#puts l.size
	
	if l.match(/^\s*>/)
		# > of fasta file as an indicator for processing last seq
		if started >0 
			puts "#{name}:#{pointer_left}:#{pointer_right}:#{seq_left}:#{width}"
			pointer_left=pointer_right+1
			pointer_right=pointer_left+l.bytesize-1
			
		else
			started =1;
			pointer_left=0;
			pointer_right=pointer_left+l.bytesize-1
			
		end
		name = l.chomp
		name = name[1,name.size]
		seq_left = pointer_right +1
		width =0
	else
		if started
			pointer_right = pointer_right + l.bytesize
			if width>0
				next
			else
				width = l.bytesize-1 # each line width of the seq
			end
		else
			next
		end
	end	
end

if started >0
	puts "#{name}:#{pointer_left}:#{pointer_right}:#{seq_left}:#{width}"
end
f.close;