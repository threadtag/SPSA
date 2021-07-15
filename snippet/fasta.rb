# Author: Hengyi Jiang <hengyi.jiang@gmail.com> <jhy@fudan.edu.cn>
# purpose: define some simple classes for fasta file parsing
# 2018.05.11
# last updated 2021-04-04
# usage: require('./fasta.rb')
# parser = Fasta_parser.new(file_name)
# seq = parser.next_fasta

require 'strscan'
require 'digest'
class Fasta_seq
	def initialize (h,s)
		@header=h
		@seq = s
	end

	def find(needle)
		r=[]
		scanner =StringScanner.new(@seq)
		while scanner.scan_until(needle) do 
			r << [scanner.pos-scanner[0].size,scanner.pos-1]
		end
		return r
	end

	def to_dna
		@seq.gsub!(/U/i,"T")
	end

	def to_rna
		@seq.gsub!(/U/i,"U")
	end

	def sanitize(type="nu")
		if (type=="nu")
			@seq.gsub!(/[^ATCGUN]/i,"")
		else
			@seq.gsub!(/[\W]/,"")
		end
	end

	attr_accessor :seq, :header
end

class Fastae_seq < Fasta_seq
	def initialize(h,s,n)
		super(h,s)
		@note=n
	end
	attr_accessor :seq, :header,:note
end


class Fasta_collection
	include Enumerable
	def initialize(file)
		@file=file
		@fasta_array = []
		@index_hash={}
	end

	def load!
		fh=Fasta_parser.new(@file)		
		i=0
		while fasta=fh.next_fasta
			md5 = Digest::MD5.hexdigest(fasta.header)
			@index_hash[md5]=i
			@fasta_array << fasta
			i +=1
		end
		fh.close
		true
	end

	def get(header)
		md5 = Digest::MD5.hexdigest(header)
		if @index_hash.has_key?(md5)
			return @fasta_array[@index_hash[md5]]
		else
			nil
		end
	end

	def add!(f)
		@fasta_array<<f
		md5 = Digest::MD5.hexdigest(f.header)
		@index_hash[md5]=@fasta_array.size-1		
	end

	def each
		@fasta_array.each do |fobj|
			yield fobj
		end
	end
	def size
		@fasta_array.size
	end
end

class Fasta_parser
	def initialize(file_name)
		begin
			@fh = File.open(file_name)
			@last_line = @fh.gets
		 	until(@last_line.match(/^\s*>/) )
		 		@last_line = @fh.gets
		 	end
		rescue
			puts ("fasta file open failed #{file_name}")
		end 
	end

	def next_fasta
		if @last_line==nil
			return false
		end
		header =@last_line.chomp
		seq =""
		while l = @fh.gets
			if(l.match(/^>/))
				@last_line=l
				break
			else
				seq << l.chomp
			end
		end
		
		if(seq !="")
			return Fasta_seq.new(header,seq)	
		else
			return false
		end						
	end

	def fasta_rewind
		@fh.rewind
		@last_line = @fh.gets
	 	until(@last_line.match(/^>/) )
	 		@last_line = @fh.gets
	 	end
	end
	def next
		next_fasta
	end
	
	def fasta_close
		@fh.close
	end
	def rewind
		fasta_rewind
	end
	def close
		@fh.close
	end
end
class Fastae_parser < Fasta_parser
	def initialize(file)
		super(file)		
	end
	def next_fastae
		if @last_line==nil
			return false
		end
		header =@last_line.chomp
		seq =""
		note=""
		while l = @fh.gets
			if(l.match(/^>/))
				@last_line=l
				break
			elsif l.match(/^#/)
				note << l.chomp
			else
				seq << l.chomp
			end
		end
		
		if(seq !="")
			return Fastae_seq.new(header,seq,note)	
		else
			return false
		end						
	end
	def next
		next_fastae
	end
	def close
		super
	end
	def rewind
		super
	end
end


def reverse_complement(seq)
	replace_map ={'A'=>"T","C"=>"G","T"=>"A","G"=>"C","a"=>'t','c'=>'g','t'=>'a','g'=>'c',"\n"=>""}
	r=seq.reverse
 	r.gsub(Regexp.union(replace_map.keys), replace_map)
end

def rev_comp_length(left, right, reversely=true,allow_wobble=true)
	# reversely
	#  3'-ATCG... <==
	#  5'-AAGC... <==
	# forword
	#  ==>  3'  GCAA...
	#  ==>  5'  CGTC...
	if left =="" or right ==""
		return 0
	end
	if !reversely
		tmp = right
		right = left
		left = tmp
	end 

	left_rvc = left.reverse
	test_len = [left.size,right.size].min	
	rr = 0
	0.upto(test_len-1) do |i|
		if is_complement?(left_rvc[i].upcase, right[i].upcase, allow_wobble)
			rr = rr +1
		else
			break
		end
	end	
	return rr
end

def is_complement?(x,y,allow_wobble=false)
	x.upcase!
	y.upcase!
	if allow_wobble
		pairs=["GC","CG","GT","TG","AT","TA","AU","UA"]
	else
		pairs=["GC","CG","AT","TA","AU","UA"]
	end
	if pairs.include?("#{x}#{y}")
		true
	else
		false
	end	
end

def is_wobble?(a,b)
	if ["GT","TG","UG","GU"].include?("#{a}#{b}".upcase)
		return true
	else
		return false
	end
end


def is_paired?(u,v,allow_wobble=true)
	# u.size == v.size assumed
	# v is in the reverse direction
	if u.size != v.size
		return false
	end
	for i in (0..(u.size-1))
		if !is_complement?(u[i],v[i],allow_wobble)
			#puts  " #{u[i]},#{v[i]}"
			return false
		end
	end
	true
end

def count_GC(s)
	c=0
	s.each_char do |ss|
		if ss=="G" or ss=="C"
			c +=1
		end
	end
	return c
end

def max_complementary(x,y)
	# x 5' ............... 3'
	#         ||||| ||
	# y 3'    ........     5'
	swapped=false
	if x.size < y.size
		# swap
		x,y =y, x
		swapped=true
	end

	min_len=4
	result=[]
	if (y.size>=min_len)
		for i in ((min_len-1)..(x.size+y.size-min_len)) do
		#        i                                                 i
		#        |                                                 |	
		#5'    ..................      ==>  5'   ..............****####
		#  ........ 5                  ==>  3'                 ........  
			case i
			when (min_len-1)..(y.size-1)
				tester_x=x[0..i]
				tester_y=y[0..i].reverse
			when (y.size-1)..(x.size-1)
				tester_x=x[(i-y.size+1)..i]
				tester_y=y.reverse
			when (x.size-1)..(x.size+y.size-min_len)
				blank_len = i+1-x.size
				overlap_len = y.size-blank_len
				tester_x=x[(x.size-overlap_len-1)..-1]
				tester_y=y[(blank_len-1)..-1].reverse
			else
				#
			end

			count=0
			continuous=0
			last_paired=-2
			max_continuous=0
			(0..(tester_x.size-1)).each do |j|
				if(is_complement?(tester_x[j],tester_y[j],true))
					#print "*"
					count +=1
					if last_paired+1 == j
						if continuous==0
							continuous=2
						else
							continuous +=1
						end
						if continuous > max_continuous
							max_continuous = continuous
						end
					else
						continuous = 0
					end
					last_paired = j
				
				else
					#print "."
				end
				
			end

			result<<[max_continuous,count,i,swapped]
		end
	else
		for i in ((y.size-1)..(x.size-1)) do
			tester_x=x[(i-y.size+1)..i]
			tester_y=y.reverse
			# puts tester_x
			# puts tester_y
			count=0
			continuous=0
			last_paired=-2
			max_continuous=0

			(0..(tester_x.size-1)).each do |j|
				if is_complement?(tester_x[j],tester_y[j],true)
					count +=1
					if last_paired+1 == j
						continuous +=1
						
						if continuous > max_continuous
							max_continuous = continuous
						end
					else
						continuous = 0
					end
					last_paired = j
				
				end				
			end
			result<<[max_continuous,count,i,swapped]	
		end	
	end

	max_continuous_record=result.max{|a,b| a[0]<=>b[0]}
	max_count_record=result.max{|a,b| a[1]<=>b[1]}

	
	# continuous sequence match is way more precedent in the alogorithm
	if max_continuous_record[0]*max_continuous_record[0] > max_count_record[1]
		return max_continuous_record
	elsif max_continuous_record[0]*max_continuous_record[0] == max_count_record[1]
		if  max_continuous_record[1]>=max_count_record[1]
			return max_continuous_record
		else
			return max_count_record
		end

	else
		return 	max_count_record
	end

	# return value  [max_continuous,count,i,swapped]
	#                i
	#                V
	# x 5' ............... 3'
	#         ||||| ||
	# y 3'    ........     5'

end