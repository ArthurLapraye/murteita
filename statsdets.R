#!/usr/bin/Rscript --vanilla -i

sfeaturefiles <- list.files(path="/home/arthur/programmes/python/stage/corpus/sfeatures/", pattern="*.txt", full.names=TRUE, recursive=FALSE)

data=new.env()

for (f in sfeaturefiles) { data[[ f ]] <- read.csv(f,header=TRUE,sep="\t") }

for (n in names(data)) { }


#lapply(files, function(x) {
 #   t <- read.csv(x, header=TRUE,sep="\t") # load file
    # apply function
   
    # write to file
    #write.table(out, "path/to/output", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#})
