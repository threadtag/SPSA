# purpose : to import index file into mysql database
# 2018-05-23 jhy @ fudan.edu.cn
# usage : ruby fasta_index_mysql_import.rb index_file db_user db_pass db_name
# database : fasta_db, table name fasta_info, structure as below

# create table fasta_info(
# 	fstid BIGINT UNSIGNED PRIMARY KEY AUTO_INCREMENT, 
# 	seq_name VARCHAR(255), 
# 	fasta_left BIGINT UNSIGNED, 
# 	fasta_right BIGINT UNSIGNED, 
# 	seq_left BIGINT UNSIGNED, 
# 	line_width INT UNSIGNED, 
# 	FULLTEXT ft_seqname(seq_name)
# )CHARACTER SET utf8 COLLATE utf8_general_ci engine=MYISAM

require 'mysql2'
require 'yaml'
if ARGV.size <2
	puts  "usage:ruby fasta_index_mysql_import.rb db_config index_file "
	abort
end
config_file = ARGV.shift
index_file = ARGV.shift

# change mariadb connection information in the db_config file
db_config=YAML.load(File.open(config_file))
db_user = db_config['db_user']
db_pass = db_config['db_password']
db_name = db_config['db_name']

begin
	fh = File.open(index_file)
rescue
	raise "file open failed:check #{index_file}"
end


client = Mysql2::Client.new(
	:host     => 'localhost',  # 主机
	:username => db_user,      # 用户名
	:password => db_pass,      # 密码
	:database => db_name,      # 数据库
	:encoding => 'utf8'        # 编码
)

fh.each do |line|
	unless (/(:\d+?){4}/.match(line) )
	    # check index file format
		raise StandardError, "Index file corrupted, please check\n";
	end
	m=/^(.+?):(\d+?):(\d+?):(\d+?):(\d+?)$/.match(line)
	seq_name = m[1]
	fasta_left = m[2]
	fasta_right = m[3]
	seq_left = m[4]
	line_width = m[5]
	if m=seq_name.match(/^(\S+)\s+\S*/)
		seq_id =m[1]
	else
		seq_id =seq_name[0..19]
	end

	insert_sql = "insert into fasta_info(seq_name,seq_id,fasta_left,fasta_right, seq_left, line_width) "
	insert_sql << "values(\"#{seq_name}\",\"#{seq_id}\",#{fasta_left},#{fasta_right},#{seq_left},#{line_width})"

	#puts insert_sql
	begin
		client.query(insert_sql)
	rescue
		puts "error inserting #{insert_sql}"
		next
	end	
end