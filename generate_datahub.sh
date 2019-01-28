# generate datahub for browser
# available: 
## 1. bg.gz for RNA-seq output
## 2. bedgraph for bedgraph file
## 3. narrowPeak for ATAC-seq peak file
## 4. bigWig for bigWig file

date
target=$1
location=$2

date_now=`date +"%m-%d-%y"`

cat <<EOF > temp.R
#!/usr/bin/env Rscript
args=commandArgs()
library(jsonlite)

target=args[6] # must be same with the file suffix
location=args[7]

files <- system(paste("ls -1"," ",location,"/","*",target, sep=""), intern=T)

datahub <- list()
for (i in 1:length(files)) {
    url=paste("http://brc.wustl.edu", files[i], sep="")
    single_file=strsplit(files[i], "/")[[1]][length(strsplit(files[i], "/")[[1]])]
    temp <- data.frame(target,url,single_file,"show","#3c67e8", 50)
    colnames(temp) <- c("type","url","name","mode","colorpositive","height")
    datahub <- append(datahub, list(temp))
}

# transfer to Json
capture.output(toJSON(datahub,pretty=T),file="temp_datahub.json")
EOF

Rscript temp.R $target $location

cat temp_datahub.json | sed '/\[/d' | sed '/\]/d' | sed 's/}/},/g' | cat <(echo "[") - <(echo "]") > dh_${date_now}_${target}.json \
    && rm temp_datahub.json \
    && mv dh_${date_now}_${target}.json $location

rm temp.R

cd $location
if [ "$target" == "narrowPeak" ]; then
    for file in `ls *$target`
	do
	bgzip $file
	tabix -p bed $file'.gz'
    done
    sed -i 's/"narrowPeak"/"bed"/g' dh_${date_now}_${target}.json
    sed -i 's/_peaks.narrowPeak"/_peaks.narrowPeak.gz"/g' dh_${date_now}_${target}.json
    sed -i 's/"name": "step3.4_peakcall_/"name": "/g' dh_${date_now}_${target}.json
elif [ "$target" == "bigWig" ]; then
    sed -i 's/"name": "step3.2_rmbl_/"name": "/g' dh_${date_now}_${target}.json
elif [ "$target" == "bedgraph" ]; then
    echo "currently no name prefix for bedgraph"
    for file in `ls *$target`
        do
        bgzip $file
        tabix -p bed $file'.gz'
    done
    sed -i 's/.bedgraph"/.bedgraph.gz"/g' dh_${date_now}_${target}.json
    sed -i 's/"type": .bedgraph.gz"/"type": "bedgraph"/g' dh_${date_now}_${target}.json #sed would recognize " as . so the type would also me edited in last step
elif [ "$target" == "bg.gz" ]; then
    echo "for RNA-seq pipe"
    sed -i 's/"bg.gz"/"bedgraph"/g' dh_${date_now}_${target}.json
    sed -i 's/"name": "step2.1_Star_/"name": "/g' dh_${date_now}_${target}.json
fi

date
echo "pipe done"



