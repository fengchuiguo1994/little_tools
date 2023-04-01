SOFT="loaddata.sh"
VERSION="1.0.0"

function usage {
    echo -e "usage: bash loadbedpe.sh -i inputfolder -g genome -f uploadfolder"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "$SOFT $VERSION"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--inputfolder INPUTFOLDER : inputfolder, the file must be [bed12, bedgraph, bedpe, bed, gene]"
    echo "   -g|--genome GENOME : Must be a genome that exists on the BASIC browser"
    echo "   -f|--folder FOLDER : Upload to the specified library"
    echo "   [-v|--version]: version"
    echo "   [-h|--help]: help"
    exit;
}

function version {
    echo -e "$SOFT version $VERSION"
    exit
}

if [ $# -lt 1 ]
then
    usage
    exit
fi

while [ -n "$1" ]
do 
    case $1 in
        -i|--inputfile) 
            input=$2
            shift 2
            ;;
        -g|--genome)
            genome=$2
            shift 2
            ;;
        -f|--folder)
            folder=$2
            shift 2
            ;;
        -h|--help)
            help
            exit 0
            ;;
        -v|--version)
            version
            exit 0
            ;;
        --) 
            shift
            break
            ;;
        -*) 
            echo "error: no such option $1."
            exit 1 
            ;;
        *) echo "error: no such option $1."
            exit 1
            ;;
        *) 
            break
            ;;
     esac
done


echo $input
echo $genome
echo $folder

ls $input/*.bed12 | while read sample
do
    echo "$sample upload..."
    samplename=basename $sample
    BED=`/opt/basic/_py/bin/python /opt/basic/console/table_util.py create ${genome} -l "${folder}" "${samplename}" | less -R | sed 's/  */\t/g' | cut -f 2 | head -3 | tail -1`
    # /opt/basic/_py/bin/python /opt/basic/console/table_util.py load -f ${BED} 1:chrom 2:start 3:end 4:name 9:color -i ${sample}
    /opt/basic/_py/bin/python /opt/basic/console/table_util.py load ${BED} 1:chrom 2:start 3:end 4:name 9:color -i ${sample}
    /opt/basic/_py/bin/python /opt/basic/console/track_util.py new ${BED} scls
done

ls $input/*.bedgraph | while read sample
do
    echo "$sample upload..."
    samplename=basename $sample
    COV=`/opt/basic/_py/bin/python /opt/basic/console/table_util.py create ${genome} -l "${folder}" "${samplename}" | less -R | sed 's/  */\t/g' | cut -f 2 | head -3 | tail -1`
    /opt/basic/_py/bin/python /opt/basic/console/track_util.py gen_cov max ${COV} ${sample}
done

ls $input/*.bedpe | while read sample
do
    echo "$sample upload..."
    samplename=basename $sample
    CLU=`/opt/basic/_py/bin/python /opt/basic/console/table_util.py create ${genome} -l "${folder}" "${samplename}" | less -R | sed 's/  */\t/g' | cut -f 2 | head -3 | tail -1`
    # /opt/basic/_py/bin/python /opt/basic/console/table_util.py load -f ${CLU} 1:chrom 2:start 3:end 4:chrom2 5:start2 6:end2 7:score -i ${sample}
    /opt/basic/_py/bin/python /opt/basic/console/table_util.py load ${CLU} 1:chrom 2:start 3:end 4:chrom2 5:start2 6:end2 7:score -i ${sample}
    /opt/basic/_py/bin/python /opt/basic/console/track_util.py new ${CLU} pcls
    /opt/basic/_py/bin/python /opt/basic/console/track_util.py new ${CLU} curv
done

ls $input/*.bed | while read sample
do  
    echo "$sample upload..."
    samplename=basename $sample
    BED=`/opt/basic/_py/bin/python /opt/basic/console/table_util.py create ${genome} -l "${folder}" "${samplename}" | less -R | sed 's/  */\t/g' | cut -f 2 | head -3 | tail -1`
    /opt/basic/_py/bin/python /opt/basic/console/table_util.py load ${BED} 1:chrom 2:start 3:end 4:name 5:score -i ${sample}
    /opt/basic/_py/bin/python /opt/basic/console/track_util.py new ${BED} scls
done

ls $input/*.gene | while read sample
do  
    echo "$sample upload..."
    samplename=basename $sample
    GENE=`/opt/basic/_py/bin/python /opt/basic/console/table_util.py create ${genome} "${samplename}" | less -R | sed 's/  */\t/g' | cut -f 2 | head -3 | tail -1`
    /opt/basic/_py/bin/python /opt/basic/console/table_util.py load_genes ${GENE} -i ${sample}
    /opt/basic/_py/bin/python /opt/basic/console/track_util.py gen_genes ${genome}
done