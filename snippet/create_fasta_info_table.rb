# purpose: to create a fasta_info table in a given msyql(Mariadb) database
# usage : ruby fasta_index_mysql_table.rb user_name password db_name
require 'mysql2'

if ARGV.size <3
	puts "usage : ruby fasta_index_mysql_table.rb user_name password db_name"
	raise ("Argument number should be 3")
end

db_user = ARGV.shift
db_pass = ARGV.shift
db_name = ARGV.shift

sql = <<EOF
create table fasta_info(
	fstid BIGINT UNSIGNED PRIMARY KEY AUTO_INCREMENT, 
	seq_name VARCHAR(255),
	seq_id VARCHAR(50), 
	fasta_left BIGINT UNSIGNED, 
	fasta_right BIGINT UNSIGNED, 
	seq_left BIGINT UNSIGNED, 
	line_width INT UNSIGNED, 
	INDEX idx_id(seq_id),
	FULLTEXT ft_seqname(seq_name)
)CHARACTER SET utf8 COLLATE utf8_general_ci engine=MYISAM
EOF


client = Mysql2::Client.new(
	:host     => 'localhost', # 主机
	:username => db_user,      # 用户名
	:password => db_pass,    # 密码
	:database => db_name,      # 数据库
	:encoding => 'utf8'       # 编码
)

begin
	client.query(sql)
rescue
	puts "error executing: #{sql}"
end
