for i in `cut -f 2  GDC_open_MAFs_manifest.txt`

do

echo $i

adress=`echo $i |cut -d'.' -f 4 `

filename=`echo $i |cut -f 2 |cut -d'.' -f 1-3,5-7 `

echo $adress $filename

wget -O "$filename" "https://gdc-api.nci.nih.gov/data/$adress"

done
