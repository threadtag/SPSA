// Author: Hengyi Jiang <hengyi.jiang@gmail.com>
package main
import(
    "fmt"
	"os"
	"bufio"
	"strings"
	// "io"
	"runtime"
)


func main(){
	if len(os.Args)<2 {
		fmt.Printf("Argument should be more than 1")
		os.Exit(1)
	}
	file_name := os.Args[1]
	f,err := os.Open(file_name)
	if err != nil{
		fmt.Println("open file failed")
		fmt.Printf("%v\n", err)
		os.Exit(1)
	}
	defer f.Close()
	runtime.GOMAXPROCS(runtime.NumCPU()-1) // MAX SPEED
	reader := bufio.NewReader(f)
	var scan_pointer int64 =0
	var seq_pointer int64 =0
	var header_name_line string
	var last_header_len int
	var line string 
	var width int
	for{
		line, err = reader.ReadString('\n')
		if err !=nil{
			break;
		}
		
		if line[0] =='>'{
			
			if scan_pointer !=0{ 
				// print the last
				fmt.Printf("%s:%d:%d:%d:%d\n",strings.TrimSuffix(header_name_line, "\n")[1:],seq_pointer, scan_pointer-1,seq_pointer+int64(last_header_len),width )
			}
			header_name_line=""
			header_name_line +=line
			seq_pointer= scan_pointer
			last_header_len =len(header_name_line)
			width=0
		}else{
			if width==0{
				width = len(line)-1
			}
		
		}
		scan_pointer += int64(len(line) )
	
	}
	fmt.Printf("%s:%d:%d:%d:%d\n",strings.TrimSuffix(header_name_line, "\n")[1:],seq_pointer, scan_pointer-1,seq_pointer+int64(last_header_len),width )
}