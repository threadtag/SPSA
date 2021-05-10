// version 2.0.1 2020-05-08
package main
import (
	"fmt"
	"os"
	"bufio"
	"strings"
	"strconv"
	"io"
	"bytes"
	"regexp"
	"flag"
	"sync"
	"runtime"
	"time"
)

type Index_line struct {
	seq_name string
	start int64
	stop int64
}

type Search_opt struct{
	string_match bool
	header_match bool
	case_sensitive bool
	both_direction bool
}

type Target_content struct{
	target_str *string
	target_regexp *regexp.Regexp
}

type Pool_seq struct {
	seq_name *string
	seq *string
} 
var file_name= flag.String("f","","fasta file path")
var index_file= flag.String("i","","index file path")
var target= flag.String("p","","pattern of the search target")
var out_file = flag.String("o","","output file, default")
var both_direction = flag.Bool("b",false,"search for both direction")
var case_sensitive = flag.Bool("c",false, "case sensitive when matching,default case insensitive, only for regular exprssion mode")

var string_match = flag.Bool("s", false, "don't use regular expression, direct string match, for speed")
var header_match = flag.Bool("h",false,"match for header only, default match sequence")

var n_threads = flag.Int("j",2,"concurrent goroutine to open")
var max_cpu = flag.Int("u",2,"max cpu cores to utilize")
var report_frequency = flag.Int("r",5000,"report frequency, default every 5000 sequences")

func main() {
	// collect the seq_idx
	usage :=`examples:
locate_subseq -j 20 -f ~/genome-db/Bacteria.fasta -i ~/genome-db/Bacteria.idx  -p "\w{3}CGGGACG\w{3}"
locate_subseq -j 20 -f ~/genome-db/Bacteria.fasta -s -p CGGGACG
-f big fasta file path
-i index file of the fasta file, when omitted search the .idx file int fasta file folder
-p pattern of the search target
-b search for both direction, default true
-h match for header only, default match sequence
-s use string match mode, no regular expression applied, faster
-c case sensitive when matching,default case insensitive, only for regular exprssion mode
-u max CPU core to use, default 2
-j concurrent goroutine to open, default 2
-r report frequency, default 5000
-o output file name
`
	flag.Parse()
	if len(os.Args) <3 {
		fmt.Println(usage)
		os.Exit(1)
	}
	if _, err := os.Stat(*file_name); err != nil {
		if os.IsNotExist(err) {
			fmt.Printf("error: fasta file:%s not found\n",*file_name)
			fmt.Println(usage)
			os.Exit(1)
		} 
	} 
	opt := Search_opt{*string_match,*header_match,*case_sensitive,*both_direction} 
	var needle Target_content
	var err error 
	if !*string_match{
		if (*case_sensitive){
			needle.target_regexp, err = regexp.Compile(*target)
		}else{			
			needle.target_regexp, err = regexp.Compile("(?i)"+ *target)
		}
		
		if err !=nil{
			fmt.Println("regular expression is not valid, please check -p option\n")
			panic(err)
		}
	}else{
		needle.target_str=  target
	}


	var index_file_path string
	if *index_file=="" {		
		t := strings.Split(*file_name,".")
		index_file_path = strings.Join(t[0:len(t)-1],".")+".idx"
	}else{
		index_file_path = *index_file
	}

	idx,err := os.Open(index_file_path)
	if err != nil{
		fmt.Printf("error: fail to open index file: %s\nuse fasta_index to generate idx file\n", index_file_path)
		fmt.Println(usage)
		os.Exit(1)
	}
	defer idx.Close()

	if *target ==""{
		fmt.Println("error: search pattern can't be empty, please use -p option to specify\n")
		fmt.Println(usage)
		os.Exit(1)
	}

	start_time := time.Now()
	// read index file
	br := bufio.NewReader(idx)
	var seq_idx []Index_line
	for{
		
		line, err := br.ReadString('\n')
		line = strings.TrimSuffix(line, "\n")
		fields := strings.Split(line,":")
		
		switch len(fields){  
		case 5:
			seq_name:=fields[0]
			start, err_ := strconv.ParseInt(fields[3], 10, 64)
			stop, err_ := strconv.ParseInt(fields[2], 10, 64)
			if err_ ==nil{
				lobj :=Index_line{seq_name: seq_name, start: start , stop: stop}
				seq_idx = append(seq_idx,lobj)				
			}
		case 6:
			seq_name:=fields[0]
			start, err_ := strconv.ParseInt(fields[4], 10, 64)
			stop, err_ := strconv.ParseInt(fields[3], 10, 64)
			if err_ ==nil{
				lobj :=Index_line{seq_name: seq_name, start: start , stop: stop}
				seq_idx = append(seq_idx,lobj)
			}
		}

		if err == io.EOF{
			break;
		}	
	}
	
	// for header search
	if *header_match {
		for i:=0; i<len(seq_idx); i++ {
			
			if !*string_match {
				if needle.target_regexp.Match([]byte(seq_idx[i].seq_name)){
					fmt.Printf("%s_%d_%d\n",seq_idx[i].seq_name,seq_idx[i].start,seq_idx[i].stop)
				}
			}else{
				if strings.Index(seq_idx[i].seq_name,*needle.target_str) != -1 {
					fmt.Printf("%s_%d_%d\n",seq_idx[i].seq_name,seq_idx[i].start,seq_idx[i].stop)
				}
			}
		}
		os.Exit(0)
	}

	// for output file
	var fout io.Writer
	if *out_file !=""{
		out, err := os.OpenFile(*out_file, os.O_WRONLY|os.O_CREATE, 0666)
		fout = out
		if err !=nil{
			fmt.Fprint(os.Stderr,"error to open output file \n")
			os.Exit(1)
		}
		defer out.Close()
	}else{
		fout=os.Stdout
	}

	// 同步开关
	// for synchronization
	var wg sync.WaitGroup 
	
	// 两个channel，一个用来放置工作项，一个用来存放处理结果。
	// two channels, one for job assigning, another for job finish signal
	ch_buffer_size := *n_threads
	runtime.GOMAXPROCS(*max_cpu)
	jobs := make(chan *Pool_seq, ch_buffer_size)
	results := make(chan int, ch_buffer_size)
	// 开启n_threads个go routine来开始计算
	// for job asignment
	for w := 1; w <= int(*n_threads); w++ {
		// func worker(id int, tgt Target_content, opt Search_opt, jobs <-chan *Pool_seq, results chan<- int)
		go worker(w, needle,  opt, jobs, results, fout, &wg)
	}

	// 开启一个 goroutine 来收集任务完成通道的信号
	// open the go routine to receive signal on result channels
	finished :=0
	to_terminate :=false
	go func(){		
		wg.Add(1)
		defer wg.Done()
		for {
			if finished >= len(seq_idx){
				//完成任务后关闭Channel
				close(jobs)
				close(results)				
				break
			}
			//  如果有返回信号, 则完成任务计数累加, 检查是否全部完成
			//  if there is finished job, increase the finished variable, check if all finished
			select{
			case <-results:
				finished ++	
			}
	
			// if finished % 5000==0{
			// 	fmt.Fprintf(os.Stderr, "[%d %%] %d sequences processed\n",100*finished/len(seq_idx),finished);
			// }
		}
		to_terminate=true
		fmt.Fprintf(os.Stderr, "[100%%] %d sequences processed\n",finished);
	}()
	

	go func(){
		// 在屏幕上监视进度
		// monitor the progress on screen
		wg.Add(1)
		defer wg.Done()
		k:=1
		last_k :=0
		for {
			if *report_frequency >0 {
				if finished >= (*report_frequency * k){
					if last_k != k{
						fmt.Fprintf(os.Stderr, "[%d %%] %d/%d  sequences processed\n",100*finished/len(seq_idx),finished,len(seq_idx) );
						last_k = k
					}
					k++
				}
			}
			if to_terminate{
				time.Sleep(2 * time.Second)
				break
			}
			// check the result every 2 second to save CPU time
			// 每2秒钟检查一下进度
			time.Sleep(2 * time.Second)						
		}	
	}()

    // 分发任务 dispatch jobs
	for i:=0; i<len(seq_idx);i++ {
		seq := read_at(*file_name, seq_idx[i].start, seq_idx[i].stop)
		jobs<- &Pool_seq{seq_name:&(seq_idx[i].seq_name),seq:seq}
	}

	fmt.Fprintf(os.Stderr,"All jobs assigned\n")
	wg.Wait()
	stop_time := time.Now()
	fmt.Fprintln(os.Stderr,"computing time:",stop_time.Sub(start_time))
	fmt.Fprintf(os.Stderr,"Controller process stopped\n")

}


//这个是工作线程，处理具体的业务逻辑，将jobs中的任务取出，处理后将处理结果放置在results中。
//This is real search process, it take our job from jobs channel, after finishing, put a signal in results channel 
func worker(id int, tgt Target_content, opt Search_opt, jobs <-chan *Pool_seq, results chan<- int, fout io.Writer,wg *sync.WaitGroup) {
	wg.Add(1)
	defer wg.Done()
	for j := range jobs {
		if !opt.string_match{
			sub_seq_regexp(fout,j.seq_name, j.seq, tgt.target_regexp, true )
			if opt.both_direction{
				sub_seq_regexp(fout,j.seq_name, j.seq, tgt.target_regexp, false)
			}
		}else{
			sub_seq_index(fout,j.seq_name, j.seq,tgt.target_str,true)
			if opt.both_direction{
				sub_seq_index(fout,j.seq_name,j.seq, tgt.target_str,false)
			}
		}
		// 清空引用,可以释放内存吗?
		// try to release the memory, will this work?
		j.seq = nil
		j.seq_name=nil
		j = nil	
		results <- id	
	}
	fmt.Fprintf(os.Stderr,"returning from worker id:%d stopped\n",id)
}



func sub_seq_regexp(fout io.Writer, header_name *string,seq *string, tgt_regexp *regexp.Regexp,forward bool){
	var pt_seq *string
	var rev_seq string
	if !forward {
		rev_seq=reverse_complement(seq)
		pt_seq = &rev_seq
	}else{
		pt_seq=seq
	}
	
	position := tgt_regexp.FindAllStringSubmatchIndex(*pt_seq,-1)
	for i:=0;i<len(position);i++{
		start_pos :=position[i][0]
		stop_pos :=position[i][1]-1
		found_seq :=(*(pt_seq))[start_pos:stop_pos+1]
		if forward{
			fmt.Fprintf(fout,">%s_fw_%d_%d\n%s\n",*header_name,start_pos+1,stop_pos+1,found_seq)
		}else{
			seq_len := len(*pt_seq)
			fmt.Fprintf(fout,">%s_rv_%d_%d\n%s\n",*header_name,seq_len-start_pos,seq_len-stop_pos,found_seq)
		}
	}
}


func sub_seq_index(fout io.Writer,header_name *string, seq *string ,t *string,forward bool) {
	var pt_seq *string
	var rev_seq string
	if !forward {
		rev_seq=reverse_complement(seq)
		pt_seq = &rev_seq
	}else{
		pt_seq=seq
	}
	pos_t := strings.Index(*pt_seq, *t)
	if pos_t == -1{
		return 
	}
	new_pos :=0
	for(pos_t != -1){
		new_pos += pos_t +len(*t) 
		if forward{
			fmt.Fprintf(fout,">%s_fw_%d_%d\n%s\n",*header_name,pos_t+1,pos_t+len(*t),*t)
		}else{
			seq_len := len(*pt_seq)
			fmt.Fprintf(fout,">%s_rv_%d_%d\n%s\n",*header_name,seq_len-pos_t,seq_len-pos_t-len(*t)+1,*t)
		}
		test_seq :=(*pt_seq)[new_pos:]
		pos_t = strings.Index(test_seq,*t)
	}	
}


func complement(s *byte) byte{
	switch {
		case *s=='C',*s=='c':
			return 'G'
		case *s=='T', *s=='t' :
			return 'A'
		case *s=='A', *s=='a' :
			return 'T'
		case *s=='G', *s=='g' :
			return 'C'
		case *s=='N', *s=='n' :
			return 'N'
		default:
			return 'X'
	}
}

func reverse_complement(s *string) string {
    r := []byte(*s)
	for i, j := 0, len(r)-1; i < len(r)/2; i, j = i+1, j-1 {
		t := complement(&r[j])
		r[j] = complement(&r[i])
		r[i] = t		
    }
    return string(r)
}



func read_at(file_path string, from int64, to int64) *string{
	file_info, err := os.Stat(file_path); 
	r :=""
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading file infomation")
		return &r
	}
	if from <0 || from > to || to >= file_info.Size(){
		fmt.Fprintf(os.Stderr, "Error: please make sure from and to don't exceed filesize, and from < t ")
		return &r
	}		 
	
	reader, error := os.Open(file_path);
    if error != nil {
        fmt.Fprintf(os.Stderr, "Error opening file %s", file_path)
	}
	defer reader.Close()
	
	var s bytes.Buffer
	var buf_size int
	if  to-from >1024{
		// each chunk is 1024 for maximum 
		buf_size =1024
	}else{
		// small span use less buffer
		buf_size=int(to-from+1)
	}
	buf := make([]byte, buf_size);
	ix := from;
	var offset int // how many bytes to read out from buf
	
	
    for {
		//ReadAt从指定的偏移量开始读取，不会改变文件偏移量
		len, ok:= reader.ReadAt(buf, int64(ix));
		if len ==0{
			break
		}
		ix = ix + int64(len);

		if (ix>=to){
			// ReadAt 每次读buf_size长，到最后一次的时候，会产生越界
			// ReadAt read buf_size bytes, at the last read for this chunk, it will cross the border
			// the whole text   ======*==========
			// the text in buf  ############
			//                        |-->to
			//                              |-->ix
			offset =int(buf_size-int((ix-to)) ) 
			if len<buf_size{
				// defend against last line EOF
				offset = len
			}
		}else{
			offset =buf_size;
		}
		// fmt.Fprintf(os.Stderr,"offset:%d\n",offset)
		b :=string(buf[0:offset])		
		s.WriteString( strings.Replace( b,  "\n","",-1)   );	// buf[offset] not included	

		if ok !=nil{
			break
		}
		if (ix>=to){
			break;
		}
	}
	r = s.String()
	return &r
}