// Author: Hengyi Jiang <hengyi.jiang@gmail.com>
package main

import (
   "fmt"
   "io/ioutil"
   "os"
   "os/exec"
   "flag"
   "strconv"
   "bufio"
   "strings"
   "io"
   "time"
   "runtime"
   "sync"
)


var rnabob_path = flag.String("b", "./rnarobo", "path of the rnabob executable bin file, default ./rnarobo")
var sub_cat_path = flag.String("s","./sub_cat","path of the sub_cat executable bin file, default./sub_cat")
var des_file = flag.String("d","","path of RNA model des file")
var file_name = flag.String("f","","path of the fasta file to be tesed")
var index_file = flag.String("i","","path of the fasta index file")
var complement_search = flag.Bool("c",false,"search both strands")
var n_threads=flag.Int("j",2,"cocurrent threads to run")
var min_len=flag.Int64("l",0,"minimum length of the each chunk")
var report_frequency = flag.Int("r",5000,"report frequency, default every 5000 chunks")

type Index_line struct {
	seq_name string
	start int64
	stop int64
}

type Search_opt struct{	
	bin string       // rnabob bin path
	sub_bin string   // suc_cat bin path
	fna  string      // big fasta file path
	des string       // rna model desc file
	complement bool  // seach both strands
}



func main() {
	
   usage :=`usage: go_rnabob -b /usr/bin/rnarobo -s /usr/bin/sub_cat -d ./trna.des -f ./bac.fasta -c -j 4
	-b path of the rnabob executable bin file, default =./rnabob
	-s path of the sub_cat executable bin file, default = ./sub_cat
	-c search both strands, default = false
	-d path of RNA model des file
	-f path of the fasta file to be tesed
	-i path of the fasta index file
	-l minimum length of the each chunk
	-j concurrent threads to run
	-r report frequency, default every 5000 chunks
note: for rnabob -F option assumed
      when the program froze, use "z + Enter" to terminate and see debug information
`
   flag.Parse()

   if len(os.Args) <2 {
		fmt.Printf(usage)
		os.Exit(1)
   }
   if _, err := os.Stat(*rnabob_path); err != nil {
		if os.IsNotExist(err) {
			fmt.Fprintf(os.Stderr,"error: bin file:%s not found\n",*rnabob_path)
			fmt.Fprintf(os.Stderr,usage)
			os.Exit(1)
		} 
	} 

   if _, err := os.Stat(*sub_cat_path); err != nil {
		if os.IsNotExist(err) {
			fmt.Fprintf(os.Stderr,"error: bin file:%s not found\n",*sub_cat_path)
			fmt.Fprintf(os.Stderr,usage)
			os.Exit(1)
		} 
	} 
	
	file_info, err := os.Stat(*file_name)
	if err != nil {
		if os.IsNotExist(err) {
			fmt.Fprintf(os.Stderr,"error: fasta file:%s not found\n",*file_name)
			fmt.Fprintf(os.Stderr,usage)
			os.Exit(1)
		} 
   }
   
   var index_file_path string
   if *index_file=="" {		
		t := strings.Split(*file_name,".")
		index_file_path = strings.Join(t[0:len(t)-1],".")+".idx"
	}else{
		index_file_path = *index_file
	}

   if _, err := os.Stat(index_file_path); err != nil {
		if os.IsNotExist(err) {
			fmt.Fprintf(os.Stderr,"error: index file:%s not found\n",index_file_path)
			fmt.Fprintf(os.Stderr,usage)
			os.Exit(1)
		} 
	} 

   if _, err := os.Stat(*des_file); err != nil {
		if os.IsNotExist(err) {
			fmt.Printf("error: model description file:%s not found\n",*des_file)
			fmt.Println(usage)
			os.Exit(1)
		} 
   }

   if *report_frequency < 0 {
		fmt.Fprintf(os.Stderr,"error:report frequency should >0, now it's set to:%d",*report_frequency)
		fmt.Println(usage)
		os.Exit(1)
   }

   if *n_threads > 2* runtime.NumCPU() {
		fmt.Fprintf(os.Stderr,"too many threads is harmful to the hard drive, please limit j to < 2*(CPU count)")
		fmt.Fprintf(os.Stderr,usage)
		os.Exit(1)
   }
 
   
	start_time := time.Now()
	// read the fna index file
	seq_idx_list, err := read_fna_idx(index_file_path)
	if err !=nil{
		fmt.Fprintf(os.Stderr,"error: fail to open index file: %s\nuse fasta_index to generate idx file\n", index_file_path)
		fmt.Println(usage)
		os.Exit(1)
	}

	// determine the chunk size
	var seq_idx []Index_line
	done_already := make(map[int]bool)
	if *min_len !=0{
		var tmp_start, tmp_stop, tmp_len int64
		s:=0
		for i:=0; i<len(seq_idx_list);i++{
			if i==0{
				tmp_start = seq_idx_list[i].start
				tmp_stop = seq_idx_list[i].stop			
			}else{
				tmp_stop = seq_idx_list[i].stop	
			}
			tmp_len = tmp_stop - tmp_start +1
			if tmp_len >= *min_len{
				seq_idx =append(seq_idx,Index_line{strconv.Itoa(s),tmp_start,tmp_stop})
				s++
				if i<len(seq_idx_list)-1{
					tmp_start=seq_idx_list[i+1].start
					tmp_stop = seq_idx_list[i+1].stop	
				}
			}else{
				if i==len(seq_idx_list)-1{
					seq_idx =append(seq_idx,Index_line{strconv.Itoa(s),tmp_start,tmp_stop})
					s++
				}
			}
			
		}
	}else{
		seq_idx =seq_idx_list
	}
	if len(seq_idx) > 0 {
		// defend against EOF
		if seq_idx[len(seq_idx)-1].stop >= file_info.Size()-1{
			seq_idx[len(seq_idx)-1].stop = file_info.Size()-2
		}

	}

	//  多个routine同步开关
	//  for synchronization of go_routines
	var wg sync.WaitGroup 
   
    // 两个channel，一个用来放置工作项，一个用来存放处理结果。
	// two channels, one for job assigning, another for job finishing signal
	const ch_buffer_size=30 
	jobs := make(chan Index_line, ch_buffer_size)
	results := make(chan int, ch_buffer_size)
 
	opt := Search_opt{	*rnabob_path,
						*sub_cat_path,
						*file_name,
						*des_file,
						*complement_search } 
	


	runtime.GOMAXPROCS(*n_threads)

	for w := 1; w <= int(*n_threads); w++ {
		go worker(w, opt ,jobs, results, &wg)

	}
	fmt.Fprintf(os.Stderr,"%d go routines running\n",*n_threads)

	// 计算第一次分发的任务 
	// 如果任务总数小于 通道缓冲数目， 第一次就全部发放出去
	// 如果任务总数大于 通道缓冲数目， 第一次将通道缓冲全部占满
	/// estimate how many routines to run in the first round
	// if job count is less than channel count, then start all the jobs in first round
	// if job count is more than channel count, send as many jobs as the channel can take
	// var first_loop_n int
	// if len(seq_idx) > ch_buffer_size{
	// 	first_loop_n = ch_buffer_size
	// }else{
	// 	first_loop_n = len(seq_idx)
	// }
	 // 分发任务, dispatching the jobs
	// for i:=0; i<first_loop_n;i++ {
	// 	jobs <- seq_idx[i]
	// }
	fmt.Fprintf(os.Stderr,"%d chunks in the waiting list\n",len(seq_idx))

  
 
	abort := make(chan struct{})
	go func(){
		z := make([]byte,1)
		for{
			os.Stdin.Read(z)
			if z[0]=='z'{
				abort <- struct{}{}
				break
			}
		}
	}()


	finished :=0
	// j :=first_loop_n // j: assigned jobs
	//获取所有的处理结果, getting the result signal
	to_terminate :=false
	j :=0 // assigned job
	go func(){
	//  开启一个线程来处理接收到的result
		wg.Add(1)
		defer wg.Done()
		for {
			if to_terminate{
				fmt.Fprintf(os.Stderr,"main process stopped, finished job count:%d\n",finished)
				// 完成任务后关闭Channel
				close(jobs)
				close(results)
				close(abort)
				break				
			}
			if finished >= len(seq_idx){
				// 完成任务后关闭Channel
				fmt.Fprintf(os.Stderr,"%d finished chunks\n",finished)
				close(jobs)
				close(results)
				close(abort)
				to_terminate=true
 				break
			}
			// 如果有线程完成计算，result在通道上有返回，继续分发任务
			// if there is a finished job, the signal in result channel will be collected 
			// and dispatch another job		
			select{
			case v,ok:=<-results:
				finished ++
				done_already[v]=true
				// fmt.Println(j,finished)
				if !ok{				
					to_terminate=true					
				}

			case <-abort:
				// for debug purpose, add an emergent abort switch
				var unfinished []int
				for ii:=0;ii<j;ii++{
					if !done_already[ii] {
						unfinished = append(unfinished,ii)
					}
				}
				fmt.Fprintf(os.Stderr, "program terminated\n")
				fmt.Fprintf(os.Stderr, "finished job count:%d, unfinished job count:%d\n",finished,len(unfinished))
				if(j<len(seq_idx)){
					fmt.Fprintf(os.Stderr, "last_job::seq_name:%s\t",seq_idx[j-1].seq_name)
					fmt.Fprintf(os.Stderr, ",start:%d\t",seq_idx[j-1].start)
					fmt.Fprintf(os.Stderr, ",stop:%d\n",seq_idx[j-1].stop)
				}
				if len(unfinished)==0{
					fmt.Fprintf(os.Stderr,"All assigned jobs finished\n")
				}else{
					fmt.Fprintf(os.Stderr,"unfinished job list:\n")
					for ii:=0;ii<len(unfinished);ii++{
						fmt.Fprintf(os.Stderr,"seq_order:%s, start:%d, stop:%d\n",seq_idx[unfinished[ii]].seq_name,seq_idx[unfinished[ii]].start,seq_idx[unfinished[ii]].stop,)
					}
				}				
				to_terminate=true							
			}
		} // end of for loop
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
				if finished > (*report_frequency * k){
					if last_k != k{
						fmt.Fprintf(os.Stderr, "%d chunks processed, total %d\n",finished,len(seq_idx));
						last_k = k
					}
					k++
				}
			}
			if to_terminate{
				time.Sleep(2 * time.Second)
				break
			}
			time.Sleep(2 * time.Second)						
		}	
	}()

	// assigning the jobs
	// 分发任务	
	for i:=0; i<len(seq_idx);i++ {
		jobs <- seq_idx[j]
		j++		
	}

	wg.Wait()
	stop_time := time.Now()
	fmt.Fprintln(os.Stderr,"computing time:",stop_time.Sub(start_time))
	fmt.Fprintf(os.Stderr,"Controller process stopped\n")
} 

func worker(id int,  opt Search_opt, jobs <-chan Index_line, results chan<- int,wg *sync.WaitGroup) {
	wg.Add(1)
	defer wg.Done()
	var j Index_line

	for  j = range jobs {
		rr,_:=strconv.Atoi(j.seq_name)
		getseq_cmd := exec.Command(opt.sub_bin,opt.fna,strconv.FormatInt(j.start,10),strconv.FormatInt(j.stop,10))
		seq, err := getseq_cmd.Output()
		if err != nil {
		   
		   fmt.Fprintf(os.Stderr,"error from sub_cat:%s, skipped\n",err.Error())  
		   // panic(err)
		   results <- rr	// make the main app coninue without loosing count
		   continue
		}
		
		if len(seq)==0 || seq[0] !='>'{
			//fmt.Fprintf(os.Stderr,"the first of the seq:%c\n", seq[0])
			fmt.Fprintf(os.Stderr,"The format of the sequence is problematic,skipped:%d~%d\n",j.start,j.stop)
			results <- rr	// make the main app coninue without loosing count
			continue
		}
		// fmt.Printf("the seq is :%s\n",seq)

		bobCmd := exec.Command(opt.bin,"-F",opt.des,"-")
		if opt.complement{
		   bobCmd = exec.Command(opt.bin,"-F","-c",opt.des,"-")
		}
		pos_t := strings.Index(opt.bin, "rnarobo")
		if pos_t != -1{
			if opt.complement{
				bobCmd = exec.Command(opt.bin,"-c",opt.des,"-")
			}else{
				bobCmd = exec.Command(opt.bin,opt.des,"-")
			}
		}
		
		bobIn, _ := bobCmd.StdinPipe()
		bobOut, _ := bobCmd.StdoutPipe()
		bobOutErr,_ :=bobCmd.StderrPipe()
		// fmt.Println("Before start")
		bobCmd.Start()
		// fmt.Println("After start")
		//bobIn.Write([]byte(">dummmy\nAAAAAAAA\n"))
		//fmt.Printf()
		bobIn.Write(seq)
		//bobIn.Write(seq[6265171:7286501+1])
		//fmt.Printf("seq-len:%d",len(seq))
		//bobIn.Write([]byte(">dummmy\nAAAAAAAA"))
		// fmt.Println("After write")

		bobIn.Close()
		bobBytes, _ := ioutil.ReadAll(bobOut)
		bobBytesErr,_:= ioutil.ReadAll(bobOutErr)
		// fmt.Println("before wait")
		bobCmd.Wait()
		// fmt.Println("After wait")
		if len(bobBytesErr) !=0{
			results <- id	// make the main app coninue without loosing count
			fmt.Fprintf(os.Stderr, "go-routine-id:%d,error-from-PIPE:%s\n",id,string(bobBytesErr))
			fmt.Fprintf(os.Stderr,"seq_name:%s,start:%d,stop:%d\n",j.seq_name,j.start,j.stop)
			fmt.Fprintf(os.Stderr,"you can use \"%s %s %d %d\" to check what sequence is causing problem\n",
			    opt.sub_bin,opt.fna ,j.start,j.stop)
			fmt.Fprintf(os.Stderr,"skipping this chunk\n")
			//panic("this goroutine shutdown\n")
			continue
		}
	
		fmt.Println(string(bobBytes))		
		results <-  rr
	}
	fmt.Fprintf(os.Stderr,"worker id:%d stopped\tj-seq_name:%s\n",id,j.seq_name)
}

func read_fna_idx(index_file_path string)([]Index_line,error){
	idx,err := os.Open(index_file_path)
	var seq_idx []Index_line
	if err != nil{
		return seq_idx,err
	}
	defer idx.Close()
   
    br := bufio.NewReader(idx)
	i :=0	
	for{
		line, err := br.ReadString('\n')
		line = strings.TrimSuffix(line, "\n")
		fields := strings.Split(line,":")
		
		switch len(fields){  
		case 5:
			// seq_name:=fields[0]
			seq_name :=strconv.Itoa(i)
			start, err_ := strconv.ParseInt(fields[1], 10, 64)
			stop, err_ := strconv.ParseInt(fields[2], 10, 64)
			if err_ ==nil{
				lobj :=Index_line{seq_name: seq_name, start: start , stop: stop}
				seq_idx = append(seq_idx,lobj)				
			}
		case 6:  //in case some fasta header contains ":"
			seq_name:=fields[0]
			start, err_ := strconv.ParseInt(fields[2], 10, 64)
			stop, err_ := strconv.ParseInt(fields[3], 10, 64)
			if err_ ==nil{
				lobj :=Index_line{seq_name: seq_name, start: start , stop: stop}
				seq_idx = append(seq_idx,lobj)
			}
		}

		if err == io.EOF{
			break;
		}
		i++	
	}
	return seq_idx,err
}