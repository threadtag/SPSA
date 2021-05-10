// Author: Hengyi Jiang <hengyi.jiang@gmail.com>
package main
import(
	"fmt"
	"os"
	"strconv"
	"bytes"
	"math/rand" 
)

	func main(){
	bases :="ATCG"
	alphabet :="ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"

	n, _err := strconv.ParseInt(os.Args[1],10,0)
	if _err !=nil{
		fmt.Println("sequence count should be number")
		os.Exit(1)
	}

	size, _err := strconv.ParseInt(os.Args[2],10,0)
	if _err !=nil{
		fmt.Println("sequence size should be number")
		os.Exit(1)
	}




	base_len :=4
	alphabet_len:=len(alphabet)

	for j:=0; j<int(n); j++{
		var buf bytes.Buffer
		var fasta_name bytes.Buffer
		for i:=0;i<20;i++{
			rand_i:=rand.Intn(alphabet_len)
			fasta_name.WriteString(alphabet[rand_i:rand_i+1])
		}		
		for i:=0;i<int(size);i++{
			rand_i:=rand.Intn(base_len)
			buf.WriteString(bases[rand_i:rand_i+1])
		}
		fmt.Printf(">%s\n%s\n",fasta_name.String(),buf.String())
	}

}