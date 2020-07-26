BEGIN {
    FS=OFS="\t"
    }

{for(i=12;i<NF;i++){
        if($i ~ "CB:Z:"){
            print $i
        }
    }
}