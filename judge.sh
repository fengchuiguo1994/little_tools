find . -type f -size +512M | while read line
do
abspath=`readlink -f $line`
filetype=`file $abspath | perl -lane 'if(/ASCII text/){print "yes"}else{print "no"}'`
# if [ -z "$typestr" ]; then
#   ls -lh $abspath
# fi
if [ "$filetype" = "yes" ];then
  ls -lh $abspath
fi
done
