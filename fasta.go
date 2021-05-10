// Author:Hengyi Jiang <hengyi.jiang@gmail.com>
// 2021-05-10
package fasta
import(
	"bufio"
	"io"
	"bytes"
	"strings"
)
type Fasta_reader struct{
	reader *bufio.Reader 
	err error
	count int  // number of seq returned until now
	last_line string
}

func New_fasta_reader(w io.Reader) Fasta_reader{
	var fr Fasta_reader
	fr.reader= bufio.NewReader(w)
	return fr
}

type Fasta_seq struct{
	Header string
	Seq string
}
type Fasta_seq_list []Fasta_seq

func (self *Fasta_reader) Next_fasta() (Fasta_seq, error){
	var buf bytes.Buffer
	if (*self).err ==nil{
		for{
			line, err := (*self).reader.ReadString('\n')
					
			line = strings.TrimSuffix(line, "\n")
			if  err == io.EOF{
				header :=(*self).last_line
				(*self).err = io.EOF				
				(*self).count +=1
				return Fasta_seq{header,buf.String()+line} , nil
			}

			if line[0]=='>'{
			
				header :=(*self).last_line
				(*self).last_line = line
				if header !=""{
					(*self).count +=1
					r := Fasta_seq{header,buf.String()}
					buf.Reset()
					return  r , nil
				}
			}
			buf.WriteString(line)
							
		}
	}	
	return Fasta_seq{"",""}, io.EOF
}





func min(x, y int) int {
    if x < y {
        return x
    }
    return y
}

func max(x,y int)int {
	if x < y {
        return y
    }
    return x
}

func complement(s byte) byte{
	switch {
		case s=='C',s=='c':
			return 'G'
		case s=='T', s=='t' :
			return 'A'
		case s=='A', s=='a' :
			return 'T'
		case s=='G', s=='g' :
			return 'C'
		case s=='N', s=='n' :
			return 'N'
		default:
			return 'X'
	}
}

func reverse(s *string) string{
	r := []byte(*s)
	for i, j := 0, len(r)-1; i < len(r)/2; i, j = i+1, j-1 {
		r[i],r[j] = r[j],r[i]		
    }
    return string(r)
}

func reverse_complement(s *string) string{
	r := []byte(*s)
	for i, j := 0, len(r)-1; i < len(r)/2; i, j = i+1, j-1 {
		r[i],r[j] = complement(r[j]),complement(r[i])		
    }
    return string(r)
}

func is_paired(a string,b string) bool{
	if len(a) !=len(b){
		return false
	}
	size := len(a)

	r:=true
	for i:=0;i<size; i++{
		if a[i] != complement(b[i]){
			r=false
			break
		} 
	}
	return r
}

func is_complement(a byte, b byte, allow_wobble bool) bool{
	if allow_wobble{
		switch(string(a)+string(b)){
		case "AT","TA","CG","GC","TG","GT":
			return true;
		default:
			return false;	
		}
	}else{
		switch(string(a)+string(b)){
		case "AT","TA","CG","GC":
			return true;
		default:
			return false;	
		}
	}
}

func pair_count(a,b string,allow_wobble bool) int{
	
	// # what count is about? 
	// # count is all the base pairs count when aligned in accordance with the consecutives
	// #        i     
	// #  --__--_______--__
	// #    ||  |||||||  ||
	// # ---++--+++++++--++----
	// #        j
	
	count := 0
	span :=min(len(a),len(b))
	for i:=0;i<span;i++{
		if is_complement(a[i],b[i],allow_wobble){
			count++
		}
	}
	
	return count
}


func calc_energy(a string ,b string) float32{
	// # Free energy calculation of modified base-pair formation in explicit solvent: A predictive model
	// #   G-C : -5.53     U-A: -4.42
	// #   U-G : -4.45     U-U: -5.82
	// # Free Energy Calculations of Watson-Crick Base Pairing in Aqueous Solution
	// #   stacking energy: (R=Purine, Y=Pyrimidine), below data is from DNA
	// #   Apc:-2.0,  TpC :-1.1
	// #   ApA : -5.7, UpU:-1.7
	var stack byte ='X'
	var energy float32=0.0
	var span int
	if len(a) < len(b){
		span = len(a)
	}else{
		span = len(b)
	}

	for i:=0;i<span;i++{
		switch strings.ToUpper( string(a[i])+string(b[i]) ){
			case "CG","GC" :
				if stack !='X'{
					switch strings.ToUpper(string(stack)+string(a[i])) {
						case "AC","TG" :
							energy += 2.0
						case "TC" :
							energy += 1.1
						case "GC","CG" :
							energy += 3.0 
						case "GG","CC" :
							energy +=5.7
						case "AG" :
							energy +=4.5
						default:
					}
				}
				stack = a[i]
				energy +=5.53
			case "AT","TA" :
					if stack !='X'{
					switch strings.ToUpper(string(stack)+string(a[i])) {
						case "TT":
							energy +=1.7
						case "AT","GT","CA","TA" :
							energy +=2.0
						case "CT" :
							energy +=1.1
						case "AA","GA":
							energy +=4.5
						default:
					}
					}
				stack = a[i]
				energy +=4.42

			case "GT","TG":
				if stack !='X'{
					switch strings.ToUpper(string(stack)+string(a[i])) {
						case "GG" :
							energy +=5.7
						case "CG" :
							energy += 3.0
						case "AG" :
							energy += 4.5
						case "TG","GT","AT" :
							energy += 2.0
						case "TT" :
							energy += 1.7
						case "CT" :
							energy +=1.1
						default:
					} 
				}
				stack = a[i]
				energy +=4.45
			default:
				stack='X'
							
		}
	}
	return energy	
}

