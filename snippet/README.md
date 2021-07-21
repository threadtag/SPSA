## To get sequences from big fasta files
### General description
Sequence pattern searching engines like `rnabob` and `infernal` can take a big fasta file as input source. These big fasta files  are usually from combined smaller sequence files, which are downloaded from public available genome databases like the `Refseq` database. When you get multiple hits from the search engines, you may want to explore the sequences around the hits, taking glimpses of the sequence context around the hits. These scripts in the folder can help to to get sequences in the big fasta file in batch. 

### How it works
For each big fasta file, an index was made and stored in a `mysql` table. At the query time, the offsets corresponding to the sequence ID will be read out from the index, and then sequence will be read by taking reference to offset values.

### Requirement
1. Enough hard disk space
2. Mysql (or Mariadb)
3. Ruby
4. Ruby gems: bio, mysql2

### Steps:
#### Install software
1. Install Mysql and Ruby according to the system manual
2. Install the packages used by the ruby scripts as below

On debian Linux and its derived distributions
```bash
sudo apt install libmysqlclient-dev
sudo apt install ruby-dev
sudo apt install make
sudo gem install rake
sudo gem install mysql2
```
On Red Hat/CentOS and other distributions using yum
```
sudo yum install mysql-devel
```
and install the gems

#### Prepare database
In the mysql client interface
```sql
>create database refseq_db
>create user  'refseq_db_user'@'localhost' identified by 'topscrete';
>grant all privileges on refseq_db.* to 'refseq_db_user'@'localhost';
>flush privileges;
```
Here, you can change database name and password to your own ones.

#### Creat table
```bash
>ruby create_fasta_info_table.rb refseq_db_user topscrete refseq_db
```
**Change the db_user and db_password and db_name**

#### Make index
```
>ruby fasta_index.rb path_to_fasta_file > path_to_fasta.idx
```

#### Import index file to mysql database
```
>ruby fasta_index_mysql_import.rb db_config.yaml path_to_fasta.idx
```
please change the db_config.yaml with the correct db_user, db_password, and db_name infomation.

#### Make query
example:
```
>ruby fasta_snippet_mysql.rb -k db_config.yaml -p /home/database/Bacteria.fna -e 0 -g 1 -l 2 -r 3 -f 4 -u 1000 -F : list.csv > target.fasta
```
The options -e -g -l -r -f are followed by a number, this number is the offset of entry_id, genome descript, start, stop, and is_forword column, this number start from 0. 
-u 1000 means: to add 1000 nt more upstream sequence to the query results.
Detailed option explanation is as below:

```
usage:ruby fasta_snippet_mysql.rb [option]  target_file.csv
      target_file.csv is the colon seperated file
      offset starting from 0
options:
-k database config file path
-p fasta_file_path
-m margin
-u upstream to include
-d downstream to include
-e entry_id in the target file
-g genome offset in the target file
-l left offset in the target file
-r right offset in the target file
-f forward? offset in the target file
-F Field seperator, default is :
```

The target file is parse from sequence search engine. Each line could be like:
```
NT_112066.1:Callithrix jacchus:3332:3358:1
```
The last field 1 is for plus strand, and 0 for minus strand
