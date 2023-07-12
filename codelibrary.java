import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.*;

// 一个工具集
public class tk {
    public static BufferedReader readFile(String fileIN) throws IOException {
        /**
        * @title: readFile
        * @description: 自动判断读文件和压缩文件
        * @author: huang xingyu
        * @date: 2021/2/26 0026 21:53
        * @param: [fileIN]
        * @return: java.io.BufferedReader
        */
        BufferedReader filehandle;
        if(fileIN.endsWith("gz") || fileIN.endsWith("GZ") || fileIN.endsWith("gzip") || fileIN.endsWith("GZIP")){
            filehandle = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileIN))));
        } else {
            filehandle = new BufferedReader(new InputStreamReader(new FileInputStream(fileIN)));
        }
        return filehandle;
    }

    public static BufferedWriter writeFile(String fileOUT) throws IOException {
        /**
        * @title: writeFile
        * @description: 自动判断写文件和压缩文件
        * @author: huang xingyu
        * @date: 2021/2/26 0026 21:54
        * @param: [fileOUT]
        * @return: java.io.BufferedWriter
        */
        BufferedWriter filehandle;
        if(fileOUT.endsWith("gz") || fileOUT.endsWith("GZ") || fileOUT.endsWith("gzip") || fileOUT.endsWith("GZIP")){
            filehandle = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(fileOUT))));
        } else {
            filehandle = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fileOUT)));
        }
        return filehandle;
    }

    public static <V> void  buildIntervalTreeWithoutRep(Map<String, IntervalTree<V>> it, String chrom, int s, int e, V o){
        /**
        * @title: buildIntervalTree
        * @description: 构建区间树，不考虑区间中是否有重叠的记录
        * @author: huang xingyu
        * @date: 2021-03-15 16:38
        * @param: [it, chrom, s, e, o]
        * @return: void
        */
        if (!it.containsKey(chrom)){
            it.put(chrom,new IntervalTree<>());
        }
        it.get(chrom).put(s, e , o);
    }

    public static <V> void  buildIntervalTree(Map<String, IntervalTree<List<V>>> it, String chrom, int s, int e, V o){
        /**
        * @title: buildIntervalTree
        * @description: 构建区间树，重叠的记录以列表的形式存储
        * @author: huang xingyu
        * @date: 2021-03-15 16:38
        * @param: [it, chrom, s, e, o]
        * @return: void
        */
        if (!it.containsKey(chrom)){
            it.put(chrom, new IntervalTree<>());
            List<V> tmplist = new ArrayList<>();
            tmplist.add(o);
            it.get(chrom).put(s, e , tmplist);
        } else {
            if (it.get(chrom).find(s,e) == null) {
                List<V> tmplist = new ArrayList<>();
                tmplist.add(o);
                it.get(chrom).put(s, e , tmplist);
            } else {
                it.get(chrom).find(s,e).getValue().add(o);
            }
        }
    }

    public static <V> Iterator<IntervalTree.Node<V>> findOverlap(Map<String, IntervalTree<V>> it, String chrom, int s, int e){
        /**
        * @title: findOverlap
        * @description: 寻找overlap的区域
        * @author: huang xingyu
        * @date: 2021-03-15 16:39
        * @param: [it, chrom, s, e]
        * @return: java.util.Iterator<htsjdk.samtools.util.IntervalTree.Node<V>>
        * @exception:
        * @throws:
        */
        return (Iterator<IntervalTree.Node<V>>) it.get(chrom).overlappers(s, e);
    }

    public static SamReader readSam(String fileIn) {
        /**
        * @title: readSam
        * @description: 读取比对的文件，bam和sam皆可
        * @author: huang xingyu
        * @date: 2021-03-15 22:11
        * @param: [fileIn]
        * @return: htsjdk.samtools.SamReader
        */
        return SamReaderFactory.makeDefault().open(new File(fileIn));
    }

    public static String getPath () {
        // discard
        // 在实际操作中发现有错误，改用单例实现
        // 获取当前程序所在的目录（不是当前运行目录）
        String realPath = tk.class.getClassLoader().getResource("").getFile();
        File file = new File(realPath);
        return file.getAbsolutePath();
    }

    public static void runCommand (String[] str) {
        // 调用cmd/shell命令
        // tk.runCommand(new String[]{"Rscript",filepath.getPath()+"/R/distance.r",(String)parser.getOption("output")+".distance.txt",(String)parser.getOption("output")+".distance.pdf"});
        try {
            Process process = Runtime.getRuntime().exec(str);
            process.waitFor();
        } catch (Exception e) {
            System.out.println(String.format("[ERROR]: failed to run : \n%s",String.join("  ",str)));
            System.out.println("You can use the cammand to get the result late\n");
            System.out.flush();
        }
    }

    public static <N> boolean setContain (Set<N> a, Set<N> b) {
        // 判断一个集合是否完全在另一个中
        Set<N> result = new TreeSet<>();
        result.addAll(a);
        result.addAll(b);
        return result.size() == a.size();
    }
}


// 一个获取当前程序所在目录的单例
// PathSingleton filepath = PathSingleton.getInstance();
// filepath.getPath()
public class PathSingleton {
    private String path;
    private static PathSingleton instance = null;
    private PathSingleton () {
        String filePath = PathSingleton.class.getProtectionDomain().getCodeSource().getLocation().getPath();
        if (filePath.endsWith(".jar")) {
            filePath = filePath.substring(0, filePath.lastIndexOf("/")+1);
        }
        File file = new File(filePath);
        path = file.getAbsolutePath();
    }
    public static PathSingleton getInstance() {
        if (instance == null) {
            instance = new PathSingleton();
        }
        return instance;
    }
    public String getPath() {
        return path;
    }
}


// 一个仍在完善中的解析参数列表的类
// CommandArgs parser = new CommandArgs("CLATK FilterBarcode - filter the mbed file", "java -cp CLATK.jar FilterBarcode -i inputfile -o outputfile [-n 1] [-a 500] [-b regionfile]", "Filter can filter out the data of some specific chromosomes or specific regions. At the same time, some records of very high complexity or very low complexity can be removed.");
// parser.addOption("-h","--help","str","",false,"Print the usage information");
// parser.addOption("-v","--version","str","",false,"Print the version","v1.0.0");
// parser.addOption("-i","--input","str","",true,"The input file");
// parser.addOption("-o","--output","str","",true,"The output file");
// parser.addOption("-n","--min","int",1,false,"The minimum complex");
// parser.parse(args);
class Parse<T> {
    private String shortpar;
    private String longpar;
    private String type;
    private T dft;
    private String info;
    private T val;
    private boolean isrequired;
    Parse (String sp, String lp, String type, T dft,boolean isrequired, String info){
        shortpar = sp;
        longpar = lp;
        this.type = type;
        this.dft = dft;
        this.isrequired = isrequired;
        this.info = info;
    }
    public void setDft(T dft) {
        this.dft = dft;
    }

    Parse(String sp, String lp, String type, T dft,boolean isrequired, String info,T val){
        shortpar = sp;
        longpar = lp;
        this.type = type;
        this.dft = dft;
        this.isrequired = isrequired;
        this.info = info;
        this.val = val;
    }
    @Override
    public String toString(){
        return String.format("%s %s %s %s",shortpar,longpar,type,info);
    }
    public String getLongpar() {
        return longpar;
    }
    public void setLongpar(String longpar) {
        this.longpar = longpar;
    }
    public String getShortpar() {
        return shortpar;
    }
    public void setShortpar(String shortpar) {
        this.shortpar = shortpar;
    }
    public void setInfo(String info) {
        this.info = info;
    }
    public String getInfo() {
        return info;
    }
    public void setType(String type) {
        this.type = type;
    }
    public String getType() {
        return type;
    }
    public T getDft() {
        return dft;
    }
    public T getVal() {
        return val;
    }
    public void setVal(T val) {
        this.val = val;
    }
    public void setIsdft(boolean isdft) {
        this.isrequired = isrequired;
    }
    public boolean getIsrequired(){
        return isrequired;
    }
}
public class CommandArgs<T> {
    private Map<String,Parse> commandargs;
    private Map<String,String> s2lpar;
    private String title;
    private String syn;
    private String des;
    private HashSet<String> required;

    public CommandArgs(){
        title = "";
        syn = "";
        des = "";
        commandargs = new LinkedHashMap<>();
        s2lpar = new LinkedHashMap<>();
        required = new HashSet<>();
    }
    public CommandArgs(String title,String syn, String des){
        this.title = title;
        this.syn = syn;
        this.des = des;
        commandargs = new LinkedHashMap<>();
        s2lpar = new LinkedHashMap<>();
        required = new HashSet<>();
    }
    public void addOption(String par1, String par2, String type, T dft, boolean isrequired, String info){
        String tmp1 = par1.trim();
        String tmp2 = par2.trim();
        String tmp3 = type.trim();
        String tmp5 = info.trim();
        if (commandargs.containsKey(tmp1)){
            error(String.format("%s parameter appears many times",tmp1));
        }
        if (commandargs.containsKey(tmp2)){
            error(String.format("%s parameter appears many times",tmp2));
        } else {
            commandargs.put(tmp2, new Parse(tmp1, tmp2, tmp3,dft,isrequired,tmp5));
            if (isrequired){
                required.add(tmp2);
            }
        }
        if (s2lpar.containsKey(tmp1)){
            error(String.format("%s parameter appears many times",tmp1));
        } else {
            s2lpar.put(tmp1,tmp2);
        }
    }
    public void addOption(String par1, String par2, String type, T dft, boolean isrequired, String info,T val){
        String tmp1 = par1.trim();
        String tmp2 = par2.trim();
        String tmp3 = type.trim();
        String tmp5 = info.trim();
        if (commandargs.containsKey(tmp1)){
            error(String.format("%s parameter appears many times",tmp1));
        }
        if (commandargs.containsKey(tmp2)){
            error(String.format("%s parameter appears many times",tmp2));
        } else {
            commandargs.put(tmp2, new Parse(tmp1, tmp2, tmp3,dft,isrequired,tmp5,val));
            if (isrequired){
                required.add(tmp2);
            }
        }
        if (s2lpar.containsKey(tmp1)){
            error(String.format("%s parameter appears many times",tmp1));
        } else {
            s2lpar.put(tmp1,tmp2);
        }
    }
    public void parse(String[] args){
        if (args.length == 0){
            usageMsg();
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--help") || args[i].equals("-h")) {
                usageMsg();
            } else if (args[i].equals("--version") || args[i].equals("-v")) {
                versionMsg();
            } else if (!args[i].startsWith("-")){
                error(String.format("Error: %s not start with -, please check the parameter.",args[i]));
            } else if (commandargs.containsKey(args[i])) {
                if ((i + 1) == args.length) {
                    error(String.format("Error: %s requires an argument.",args[i]));
                }
                required.remove(args[i]);
                if (commandargs.get(args[i]).getType().equals("int")) {
                    try {
                        commandargs.get(args[i]).setVal(Integer.parseInt(args[++i]));
                    } catch (NumberFormatException e){
                        e.printStackTrace();
                        error(String.format("Error: %s'%s must be int.",args[i-1],args[i]));
                    }
                } else if (commandargs.get(args[i]).getType().equals("float")) {
                    try {
                        commandargs.get(args[i]).setVal(Float.parseFloat(args[++i]));
                    } catch (NumberFormatException e){
                        e.printStackTrace();
                        error(String.format("Error: %s'%s must be float.",args[i-1],args[i]));
                    }
                } else {
                    commandargs.get(args[i]).setVal(args[++i]);
                }
            } else {
                if (s2lpar.containsKey(args[i])){
                    if (commandargs.containsKey(s2lpar.get(args[i]))){
                        if ((i + 1) == args.length) {
                            error(String.format("Error: %s requires an argument.",args[i]));
                        }
                        required.remove(s2lpar.get(args[i]));
                        if (commandargs.get(s2lpar.get(args[i])).getType().equals("int")) {
                            try {
                                commandargs.get(s2lpar.get(args[i])).setVal(Integer.parseInt(args[++i]));
                            } catch (NumberFormatException e){
                                e.printStackTrace();
                                error(String.format("Error: %s'%s must be int.",args[i-1],args[i]));
                            }
                        } else if (commandargs.get(s2lpar.get(args[i])).getType().equals("float")) {
                            // System.out.println("float"+commandargs.get(s2lpar.get(args[i])).getType()+commandargs.get(s2lpar.get(args[i])));
                            try {
                                commandargs.get(s2lpar.get(args[i])).setVal(Float.parseFloat(args[++i]));
                            } catch (NumberFormatException e){
                                e.printStackTrace();
                                error(String.format("Error: %s'%s must be float.",args[i-1],args[i]));
                            }
                        } else {
                            // System.out.println("else"+commandargs.get(s2lpar.get(args[i])).getType()+commandargs.get(s2lpar.get(args[i])));
                            commandargs.get(s2lpar.get(args[i])).setVal(args[++i]);
                        }
                    } else {
                        error(String.format("Error: argument not recognized: %s",args[i]));
                    }
                } else {
                    error(String.format("Error: argument not recognized: %s", args[i]));
                }
            }
        }
        if (required.size() > 0){
            error(String.format("Error: %s argument must be recognized: ", required));
        }
    }
    public T getOption(String str){
        String tmp = "--"+str;
        if (!commandargs.containsKey(tmp)){
            throw new IllegalArgumentException(String.format("%s does not exist",tmp));
        }
        if (commandargs.get(tmp).getVal() != null){
            return (T) commandargs.get(tmp).getVal();
        } else {
            return (T) commandargs.get(tmp).getDft();
        }
    }
    @Override
    public String toString(){
        StringBuffer sbuff = new StringBuffer();
        for(String x:commandargs.keySet()){
            sbuff.append(String.format("%s: %s\n",x,commandargs.get(x).toString()));
        }
        return sbuff.toString();
    }
    public void error(String errmsg) {
        System.err.println(errmsg);
        System.exit(1);
    }
    public void usageMsg(){
        System.out.println(String.format("\t\t\t\t%s\t\t",title));
        System.out.println("SYNOPSIS");
        System.out.println(String.format("\t\t\t%s\t",syn));
        System.out.println("DESCRIPTION");
        System.out.println(String.format("\t\t%s\t",des));
        System.out.println(String.format("\tThe options for the program as follows:"));
        for(String x: commandargs.keySet()){
            if (x.equals("--help")){
                System.out.println(String.format("\t%s  %s\t\t %s ",commandargs.get(x).getShortpar(),commandargs.get(x).getLongpar(),commandargs.get(x).getInfo()));
            } else if (x.equals("--version")){
                System.out.println(String.format("\t%s  %s\t\t %s ",commandargs.get(x).getShortpar(),commandargs.get(x).getLongpar(),commandargs.get(x).getInfo()));
            } else if (!commandargs.get(x).getIsrequired()){
                System.out.println(String.format("\t%s  %s\t\t %s [%s][default: %s]",commandargs.get(x).getShortpar(),commandargs.get(x).getLongpar(),commandargs.get(x).getInfo(),commandargs.get(x).getType(),commandargs.get(x).getDft()));
            } else {
                System.out.println(String.format("\t%s  %s\t\t %s [%s]",commandargs.get(x).getShortpar(),commandargs.get(x).getLongpar(),commandargs.get(x).getInfo(),commandargs.get(x).getType()));
            }
        }
        System.exit(0);
    }
    public void versionMsg(){
        System.out.println(commandargs.get("--version").getVal());
        System.exit(0);
    }
    public void setSyn(String syn) {
        this.syn = syn;
    }
    public String getSyn() {
        return syn;
    }
    public void setTitle(String title) {
        this.title = title;
    }
    public String getTitle() {
        return title;
    }
    public String getDes() {
        return des;
    }
    public void setDes(String des) {
        this.des = des;
    }
}


// 针对生物信息学中常用文件格式的一些类（这里的region是从1开始的，为了和htsjdk的一致，所以区间是左闭右闭）
public class Region implements Cloneable {
    private int start;
    private int end;
    public Region(int s, int e){
        start = s;
        end = e;
        if (s > e){ // 序列从1开始，所以是存在相等的情况的
            throw new IllegalArgumentException("beginning of range must be less and equal than end");
        }
    }
    @Override
    public String toString(){
        return String.format("%d\t%d",start,end);
    }
    @Override
    public Object clone() throws CloneNotSupportedException {
        Region obj = null;
        try {
            obj = (Region) super.clone();
        } catch (CloneNotSupportedException ex){
            ex.printStackTrace();
        }
        return obj;
    }
    public boolean contains(int pos) { // [] 左闭右闭的区间
        return  pos <= end && pos >= start;
    }
    public boolean contains(int s, int e){
        return s <= end && s >= start && e <= end && e >= start;
    }
    public boolean contains (Region o) {
        return o.getStart() <= end && o.getStart() >= start && o.getEnd() <= end && o.getEnd() >= start;
    }
    public boolean intersect(int s, int e){
        return e >= start && s <= end;
    }
    public boolean intersect(Region o){
        return o.getEnd() >= start && o.getStart() <= end;
    }
    public int len () {
        return end - start + 1;
    }
    public int getStart() {
        return start;
    }
    public void setStart(int start) {
        this.start = start;
    }
    public int getEnd() {
        return end;
    }
    public void setEnd(int end) {
        this.end = end;
    }
}
public class Bed<V> extends Region implements Comparable<Bed>, Cloneable {
    private String chrom;
    private int start;
    private int end;
    private String name;
    private float cov;
    private V value; // 方便添加一些其它信息，所以就只有默认了bed3格式，其余的信息可以以字典或者其它类存储，然后赋值该变量
    public Bed(String c, int s, int e){
        super(s,e);
        start = s;
        end = e;
        chrom = c;
        name = null;
        cov = 0;
        value = null;
    }
    public Bed(String c, int s, int e, String n){
        super(s,e);
        start = s;
        end = e;
        chrom = c;
        name = n;
        cov = 0;
        value = null;
    }
    @Override
    public int compareTo(Bed o){
        if (this == o){
            return 0;
        }
        int result = chrom.compareTo(o.getChrom());
        if (result == 0){
            result = start - o.getStart();
            if (result == 0){
                result = end - o.getEnd();
            }
        }
        return result;
    }
    @Override
    public Object clone() throws CloneNotSupportedException{
        Bed obj = null;
        try {
            obj = (Bed) super.clone();
        } catch (CloneNotSupportedException ex){
            ex.printStackTrace();
        }
        return obj;
    }
    @Override
    public String toString(){
        return String.format("%s\t%d\t%d",chrom,start,end);
    }
    public String toString(int a) { // 为了方便输出为标准bed格式
        return String.format("%s\t%d\t%d",chrom,start-1,end);
    }
    public boolean contains(String c, int pos){
        if (chrom.compareTo(c) == 0){
            return  pos <= end && pos >= start;
        }
        return false;
    }
    public boolean contains(String c, int s, int e){
        if (chrom.compareTo(c) == 0){
            return s <= end && s >= start && e <= end && e >= start;
        }
        return false;
    }
    public boolean contains(Bed o){
        if (chrom.compareTo(o.getChrom()) == 0){
            return contains(o.getChrom(), o.getStart(), o.getEnd());
        }
        return false;
    }
    public boolean intersect(String c, int s, int e){
        if (chrom.compareTo(c) == 0) {
            return e >= start && s <= end;
        }
        return false;
    }
    public boolean intersect(Bed o){
        if (chrom.compareTo(o.getChrom()) == 0) {
            return intersect(o.getChrom(),o.getStart(), o.getEnd());
        }
        return false;
    }
    public void mergeRegion (String c, int s, int e){
        if (start > s){
            setStart(s);
        }
        if (end < e){
            setEnd(e);
        }
    }
    public int overlapDis (String c, int s, int e) {
        if (chrom.equals(c)) {
            if (intersect(c,s,e)) {
                return 0;
            } else {
                if (s > end) {
                    return s-end;
                } else {
                    return start-e;
                }
            }
        } else {
            return -1;
        }
    }
    public int overlapDis (Bed o){
        return overlapDis(o.getChrom(),o.getStart(),o.getEnd());
    }
    public void mergeRegion (Bed o) {
        mergeRegion(o.getChrom(), o.getStart(), o.getEnd());
    }
    public String getChrom() {
        return chrom;
    }
    public void setChrom(String chrom) {
        this.chrom = chrom;
    }
    public String getName() {
        return name;
    }
    public void setName(String name) {
        this.name = name;
    }
    public float getCov() {
        return cov;
    }
    public void setCov(float cov) {
        this.cov = cov;
    }
    public V getValue() {
        return value;
    }
    public void setValue(V value) {
        this.value = value;
    }
    @Override
    public int getEnd() {
        return end;
    }
    @Override
    public void setEnd(int end) {
        this.end = end;
    }
    @Override
    public int getStart() {
        return start;
    }
    @Override
    public void setStart(int start) {
        this.start = start;
    }
    @Override
    public int len() {
        return end - start + 1;
    }
}
public class Bedpe<V2> implements Comparable<Bedpe<V2>>,Cloneable{
    private Bed<V> bed1;
    private Bed<V> bed2;
    private float interaction;
    private int mark; // 该参数是用来标识的。可省略
    private V2 bedpeValue; // 跟bed一样，方便用户自定义其它的属性

    public Bedpe(Bed<V> bed1, Bed<V> bed2) {
        this.bed1 = bed1;
        this.bed2 = bed2;
        interaction = 1;
        mark = 0;
        bedpeValue = null;
    }
    public Bedpe(Bed<V> bed1, Bed<V> bed2, float interaction){
        this.bed1 = bed1;
        this.bed2 = bed2;
        this.interaction = interaction;
        mark = 0;
        bedpeValue = null;
    }
    @Override
    public String toString(){
        return String.format("%s\t%s",bed1.toString(),bed2.toString());
    }
    public String toString(int a){ // 方便输出标准格式的bed文件
        return String.format("%s\t%s",bed1.toString(1),bed2.toString(1));
    }
    @Override
    public Object clone() throws CloneNotSupportedException {
        Bedpe<V2> u = null;
        try{
            u = (Bedpe<V2>) super.clone();
            u.bed1 = (Bed<V>) bed1.clone();
            u.bed2 = (Bed<V>) bed2.clone();
        } catch (CloneNotSupportedException ex){
            ex.printStackTrace();
        }
        return u;
    }
    @Override
    public int compareTo(Bedpe<V2> o){ // 方便其中一种的排序方式
        int result = this.bed1.getChrom().compareTo((o.bed1.getChrom()));
        if(result == 0){
            result = this.bed2.getChrom().compareTo((o.bed2.getChrom()));
            if(result == 0) {
                result = this.bed1.getStart() - o.bed1.getStart();
                if (result == 0) {
                    result = this.bed1.getEnd() - o.bed1.getEnd();
                    if (result == 0) {
                        result = this.bed2.getStart() - o.bed2.getStart();
                        if(result == 0){
                            result = this.bed2.getEnd() - o.bed2.getEnd();
                        }
                    }
                }
            }
        }
        return result;
    }
    public void reverse(){
        Bed<V> tmp = this.bed1;
        this.bed1 = this.bed2;
        this.bed2 = bed1;
    }
    public Boolean sorted(){
        if (this.bed1.compareTo(this.bed2) > 0) { // big small
            this.reverse();
            return true;
        }
        return false;
    }
    public boolean intersect(Bed<V> bed1,Bed<V> bed2){
        return this.bed1.intersect(bed1) && this.bed2.intersect(bed2)
    }
    public boolean intersect(Bedpe<V2> o){
        return intersect(o.getBed1(),o.getBed2());
    }
    public boolean contains(String c1, int s1, int e1, String c2, int s2, int e2){
        if (bed1.getChrom().compareTo(c1) == 0 && bed2.getChrom().compareTo(c2) == 0){
            return s1 <= bed1.getEnd() && s1 >= bed1.getStart() && e1 <= bed1.getEnd() && e1 >= bed1.getStart() &&
                    s2 <= bed2.getEnd() && s2 >= bed2.getStart() && e2 <= bed2.getEnd() && e2 >= bed2.getStart() ;
        }
        return false;
    }
    public boolean contains(Bedpe<V2> o){
        if (bed1.getChrom().compareTo(o.bed1.getChrom()) == 0 && bed2.getChrom().compareTo(o.bed2.getChrom()) == 0){
            return contains(o.bed1.getChrom(), o.bed1.getStart(), o.bed1.getEnd(),o.bed2.getChrom(), o.bed2.getStart(), o.bed2.getEnd());
        }
        return false;
    }
    public void mergeRegion(Bed<V> bed1,Bed<V> bed2){
        this.bed1.mergeRegion(bed1);
        this.bed2.mergeRegion(bed2);
        interaction += 1;
    }
    public void mergeRegion(Bedpe<V2> o){
        this.bed1.mergeRegion(o.getBed1());
        this.bed2.mergeRegion(o.getBed2());
        this.interaction += o.getInteraction();
    }
    public int getMark() {
        return mark;
    }
    public void setMark(int mark) {
        this.mark = mark;
    }
    public float getInteraction() {
        return interaction;
    }
    public void setInteraction(float interaction) {
        this.interaction = interaction;
    }
    public Bed<V> getBed1() {
        return bed1;
    }
    public void setBed1(Bed<V> bed1) {
        this.bed1 = bed1;
    }
    public Bed<V> getBed2() {
        return bed2;
    }
    public void setBed2(Bed<V> bed2) {
        this.bed2 = bed2;
    }
    public V2 getBedpeValue() {
        return bedpeValue;
    }
    public void setBedpeValue(V2 bedpeValue) {
        this.bedpeValue = bedpeValue;
    }
}



// 多路归并排序bed和bedpe文件
// MultiMergeSort aa = new MultiMergeSort("infile");
// aa.splitBed("infile"+".tmp",20000000);
// aa.multiMergeBed("outfile");
// MultiMergeSort aa = new MultiMergeSort("infile");
// aa.splitBedpe(in+".tmp",20000000);
// aa.multiMergeBedpe("outfile");
public class MultiMergeSort {
    private String fileName;
    private List<String> outFiles = new ArrayList<>();
    public MultiMergeSort(String f) {
        fileName = f;
    }
    public void splitBed (String prefix, int n) throws IOException {
        BufferedReader fin = tk.readFile(fileName);
        int Nfile = 0;
        String line;
        int flag = 0;
        BufferedWriter fout = tk.writeFile(prefix+"."+flag+".tmp");
        outFiles.add(prefix+"."+flag+".tmp");
        String[] tmp;
        Bed bed;
        List<Bed> beds = new ArrayList<>(n);
        while ((line = fin.readLine()) != null ){
            if (line.length() <= 0){
                continue;
            }
            Nfile += 1;
            tmp = line.split("\\s+",4);
            bed = new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2]));
            if (tmp.length == 4) {
                bed.setValue(tmp[3]);
            }
            beds.add(bed);
            if (Nfile%n == 0) {
                Collections.sort(beds);
                for (Bed x: beds) {
                    fout.write(x.toString(1));
                    if (x.getValue()==null){
                        fout.write("\n");
                    } else {
                        fout.write("\t"+x.getValue()+"\n");
                    }
                }
                beds.clear();
                fout.close();
                flag += 1;
                fout = tk.writeFile(prefix+"."+flag+".tmp");
                outFiles.add(prefix+"."+flag+".tmp");
            }
        }
        if (beds.size() != 0) {
            Collections.sort(beds);
            for (Bed x: beds) {
                fout.write(x.toString(1));
                if (x.getValue()==null){
                    fout.write("\n");
                } else {
                    fout.write("\t"+x.getValue()+"\n");
                }
            }
            fout.close();
        } else {
            fout.close();
            File tmpfile = new File(prefix+"."+flag+".tmp");
            tmpfile.delete();
        }
    }

    public void multiMergeBed (BufferedWriter fout) throws IOException {
        int fileSize = outFiles.size();
        if (fileSize == 1) {
            BufferedReader fin = tk.readFile(outFiles.get(0));
            String line;
            while((line = fin.readLine()) != null){
                fout.write(line+"\n");
            }
            File tmpfile = new File(outFiles.get(0));
            tmpfile.delete();
            return;
        }
        List<BufferedReader> fins = new ArrayList<>();
        Bed[] fileCompare = new Bed[fileSize];
        for (int i = 0; i < fileSize; i++) {
            fins.add(tk.readFile(outFiles.get(i)));
        }
        int index = 0;
        String[] tmp;
        Bed bed;
        for (int i=0; i < fileSize; i++) {
            tmp = fins.get(i).readLine().split("\\s+",4);
            bed = new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2]));
            if (tmp.length == 4) {
                bed.setValue(tmp[3]);
            }
            fileCompare[i] = bed;
        }
        int count = fileSize;
        int[] sum = new int[fileSize];
        while (count > 1){
            index = getMinIndex(fileCompare);
            fout.write(fileCompare[index].toString(1));
            if (fileCompare[index].getValue()==null){
                fout.write("\n");
            } else {
                fout.write("\t"+fileCompare[index].getValue()+"\n");
            }
            sum[index]++;
            try {
                tmp = fins.get(index).readLine().split("\\s+",4);

                bed = new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2]));
                if (tmp.length == 4) {
                    bed.setValue(tmp[3]);
                }
                fileCompare[index] = bed;
            } catch (Exception e) {
                fileCompare[index] = null;
                count--;
                fins.get(index).close();
            }
        }
        String line;
        int sIndex = getSindex(fileCompare);
        fout.write(fileCompare[sIndex].toString(1));
        if (fileCompare[sIndex].getValue()==null){
            fout.write("\n");
        } else {
            fout.write("\t"+fileCompare[sIndex].getValue()+"\n");
        }
        while ((line = fins.get(sIndex).readLine()) != null) {
            fout.write(line+"\n");
        }
        fins.get(sIndex).close();
        for (String x: outFiles){
            File tmpfile = new File(x);
            tmpfile.delete();
        }
    }

    public void multiMergeBed (String output) throws IOException {
        int fileSize = outFiles.size();
        if (fileSize == 1) {
            File tmp = new File(outFiles.get(0));
            tmp.renameTo(new File(output));
            return;
        }
        BufferedWriter fout = tk.writeFile(output);
        List<BufferedReader> fins = new ArrayList<>();
        Bed[] fileCompare = new Bed[fileSize];
        for (int i = 0; i < fileSize; i++) {
            fins.add(tk.readFile(outFiles.get(i)));
        }
        int index = 0;
        String[] tmp;
        Bed bed;
        for (int i=0; i < fileSize; i++) {
            tmp = fins.get(i).readLine().split("\\s+",4);
            bed = new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2]));
            if (tmp.length == 4) {
                bed.setValue(tmp[3]);
            }
            fileCompare[i] = bed;
        }
        int count = fileSize;
        int[] sum = new int[fileSize];
        while (count > 1){
            index = getMinIndex(fileCompare);
            fout.write(fileCompare[index].toString(1));
            if (fileCompare[index].getValue()==null){
                fout.write("\n");
            } else {
                fout.write("\t"+fileCompare[index].getValue()+"\n");
            }
            sum[index]++;
            try {
                tmp = fins.get(index).readLine().split("\\s+",4);

                bed = new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2]));
                if (tmp.length == 4) {
                    bed.setValue(tmp[3]);
                }
                fileCompare[index] = bed;
            } catch (Exception e) {
                fileCompare[index] = null;
                count--;
                fins.get(index).close();
            }
        }
        String line;
        int sIndex = getSindex(fileCompare);
        fout.write(fileCompare[sIndex].toString(1));
        if (fileCompare[sIndex].getValue()==null){
            fout.write("\n");
        } else {
            fout.write("\t"+fileCompare[sIndex].getValue()+"\n");
        }
        while ((line = fins.get(sIndex).readLine()) != null) {
            fout.write(line+"\n");
        }
        fins.get(sIndex).close();
        fout.close();
        for (String x: outFiles){
            File tmpfile = new File(x);
            tmpfile.delete();
        }
    }

    // 最后一个文件
    public static <V> int getSindex(V[] x) {
        int result = 0;
        for (int i = 0; i < x.length; i++) {
            if (x[i] != null) {
                result = i;
                break;
            }
        }
        return result;
    }

    // 最小的数据
    public static <V extends Comparable> int getMinIndex(V[] x){
        int flag = -1;
        for (int i = 0; i < x.length; i++) {
            if (x[i] == null) {
                continue;
            } else if (x[i] != null && flag == -1) {
                flag = i;
            } else if (x[i] != null && flag != -1) {
                if (x[i].compareTo(x[flag]) < 0) {
                    flag = i;
                }
            }
        }
        return flag;
    }

    public void splitBedpe (String prefix, int n) throws IOException {
        BufferedReader fin = tk.readFile(fileName);
        int Nfile = 0;
        String line;
        int flag = 0;
        BufferedWriter fout = tk.writeFile(prefix+"."+flag+".tmp");
        outFiles.add(prefix+"."+flag+".tmp");
        String[] tmp;
        Bedpe bedpe;
        List<Bedpe> bedpes = new ArrayList<>(n);
        while ((line = fin.readLine()) != null ){
            if (line.length() <= 0){
                continue;
            }
            Nfile += 1;
            tmp = line.split("\\s+",7);
            bedpe = new Bedpe(new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2])),new Bed(tmp[3],Integer.parseInt(tmp[4])+1,Integer.parseInt(tmp[5])));
            if (tmp.length == 7){
                bedpe.setBedpeValue(tmp[6]);
            }
            bedpes.add(bedpe);
            if (Nfile%n == 0) {
                Collections.sort(bedpes);
                for (Bedpe x: bedpes) {
                    fout.write(x.toString(1));
                    if (x.getBedpeValue() == null) {
                        fout.write("\n");
                    } else {
                        fout.write("\t"+x.getBedpeValue()+"\n");
                    }
                }
                bedpes.clear();
                fout.close();
                flag += 1;
                fout = tk.writeFile(prefix+"."+flag+".tmp");
                outFiles.add(prefix+"."+flag+".tmp");
            }
        }
        if (bedpes.size() != 0) {
            Collections.sort(bedpes);
            for (Bedpe x: bedpes) {
                fout.write(x.toString(1));
                if (x.getBedpeValue() == null) {
                    fout.write("\n");
                } else {
                    fout.write("\t"+x.getBedpeValue()+"\n");
                }
            }
            fout.close();
        } else {
            fout.close();
            File tmpfile = new File(prefix+"."+flag+".tmp");
            tmpfile.delete();
        }
    }

    public void multiMergeBedpe (String output) throws IOException {
        int fileSize = outFiles.size();
        if (fileSize == 1) {
            File tmp = new File(outFiles.get(0));
            tmp.renameTo(new File(output));
            return;
        }
        BufferedWriter fout = tk.writeFile(output);
        List<BufferedReader> fins = new ArrayList<>(); // 文件句柄
        Bedpe[] fileCompare = new Bedpe[fileSize]; // 待比较的记录
        for (int i = 0; i < fileSize; i++) {
            fins.add(tk.readFile(outFiles.get(i)));
        }
        int index = 0;
        String[] tmp;
        Bedpe bedpe;
        for (int i=0; i < fileSize; i++) {
            tmp = fins.get(i).readLine().split("\\s+",7);
            bedpe = new Bedpe(new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2])),new Bed(tmp[3],Integer.parseInt(tmp[4])+1,Integer.parseInt(tmp[5])));
            if (tmp.length == 7){
                bedpe.setBedpeValue(tmp[6]);
            }
            fileCompare[i] = bedpe;
        }
        int count = fileSize;
        int[] sum = new int[fileSize];
        while (count > 1){  // 开始迭代输出每次最小的记录
            index = getMinIndex(fileCompare);
            fout.write(fileCompare[index].toString(1));
            if (fileCompare[index].getBedpeValue() == null) {
                fout.write("\n");
            } else {
                fout.write("\t"+fileCompare[index].getBedpeValue()+"\n");
            }
            sum[index]++;
            try {
                tmp = fins.get(index).readLine().split("\\s+",7);
                bedpe = new Bedpe(new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2])),new Bed(tmp[3],Integer.parseInt(tmp[4])+1,Integer.parseInt(tmp[5])));
                if (tmp.length == 7){
                    bedpe.setBedpeValue(tmp[6]);
                }
                fileCompare[index] = bedpe;
            } catch (Exception e) {
                fileCompare[index] = null;
                count--;
                fins.get(index).close();
            }
        }
        String line;
        int sIndex = getSindex(fileCompare);
        fout.write(String.format("%s",fileCompare[sIndex].toString(1)));
        if (fileCompare[sIndex].getBedpeValue() == null) {
            fout.write("\n");
        } else {
            fout.write("\t"+fileCompare[sIndex].getBedpeValue()+"\n");
        }
        while ((line = fins.get(sIndex).readLine()) != null) {
            fout.write(line+"\n");
        }
        fins.get(sIndex).close();
        fout.close();
        for (String x: outFiles){
            File tmpfile = new File(x);
            tmpfile.delete();
        }
    }

    public void multiMergeBedpe (BufferedWriter fout) throws IOException {
        int fileSize = outFiles.size();
        if (fileSize == 1) {
            BufferedReader fin = tk.readFile(outFiles.get(0));
            String line;
            while((line = fin.readLine()) != null){
                fout.write(line+"\n");
            }
            File tmpfile = new File(outFiles.get(0));
            tmpfile.delete();
            return;
        }
        List<BufferedReader> fins = new ArrayList<>(); // 文件句柄
        Bedpe[] fileCompare = new Bedpe[fileSize]; // 待比较的记录
        for (int i = 0; i < fileSize; i++) {
            fins.add(tk.readFile(outFiles.get(i)));
        }
        int index = 0;
        String[] tmp;
        Bedpe bedpe;
        for (int i=0; i < fileSize; i++) {
            tmp = fins.get(i).readLine().split("\\s+",7);
            bedpe = new Bedpe(new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2])),new Bed(tmp[3],Integer.parseInt(tmp[4])+1,Integer.parseInt(tmp[5])));
            if (tmp.length == 7){
                bedpe.setBedpeValue(tmp[6]);
            }
            fileCompare[i] = bedpe;
        }
        int count = fileSize;
        int[] sum = new int[fileSize];
        while (count > 1){  // 开始迭代输出每次最小的记录
            index = getMinIndex(fileCompare);
            fout.write(fileCompare[index].toString(1));
            if (fileCompare[index].getBedpeValue() == null) {
                fout.write("\n");
            } else {
                fout.write("\t"+fileCompare[index].getBedpeValue()+"\n");
            }
            sum[index]++;
            try {
                tmp = fins.get(index).readLine().split("\\s+",7);
                bedpe = new Bedpe(new Bed(tmp[0],Integer.parseInt(tmp[1])+1,Integer.parseInt(tmp[2])),new Bed(tmp[3],Integer.parseInt(tmp[4])+1,Integer.parseInt(tmp[5])));
                if (tmp.length == 7){
                    bedpe.setBedpeValue(tmp[6]);
                }
                fileCompare[index] = bedpe;
            } catch (Exception e) {
                fileCompare[index] = null;
                count--;
                fins.get(index).close();
            }
        }
        String line;
        int sIndex = getSindex(fileCompare);
        fout.write(String.format("%s",fileCompare[sIndex].toString(1)));
        if (fileCompare[sIndex].getBedpeValue() == null) {
            fout.write("\n");
        } else {
            fout.write("\t"+fileCompare[sIndex].getBedpeValue()+"\n");
        }
        while ((line = fins.get(sIndex).readLine()) != null) {
            fout.write(line+"\n");
        }
        fins.get(sIndex).close();
        for (String x: outFiles){
            File tmpfile = new File(x);
            tmpfile.delete();
        }
    }
}



public class Bed<N extends Number & Comparable,V> implements Comparable<Bed>,Cloneable{
    private String chrom;
    private N start;
    private N end;
    private String name;
    private float cov;
    private V value;
    public Bed(String c,N s, N e){
        chrom = c;
        start = s;
        end = e;
        if(start.compareTo(end) > 0 ){ // 序列从0开始
            throw new IllegalArgumentException("beginning of range must be less and equal than end");
        }
        name = null;
        cov = 0f;
        value = null;
    }
    public Bed(String c,N s, N e, String n){
        this(c,s,e);
        name = n;
    }
    public Bed(String c,N s, N e, V v){
        this(c,s,e);
        value = v;
    }
    public void setChrom(String chrom) {
        this.chrom = chrom;
    }
    public String getChrom() {
        return chrom;
    }
    public void setStart(N start) {
        this.start = start;
    }
    public N getStart() {
        return start;
    }
    public void setEnd(N end) {
        this.end = end;
    }
    public N getEnd() {
        return end;
    }
    public void setName(String name) {
        this.name = name;
    }
    public String getName() {
        return name;
    }
    public void setCov(float cov) {
        this.cov = cov;
    }
    public float getCov() {
        return cov;
    }
    public void setValue(V value) {
        this.value = value;
    }
    public V getValue() {
        return value;
    }
    @Override
    public int compareTo(Bed o) {
        int res = 0;
        if (this != o) {
            res = chrom.compareTo(o.getChrom());
            if (res != 0) {
                res = start.compareTo(o.getStart());
                if (res != 0) {
                    res = end.compareTo(o.getEnd());
                }
            }
        }
        return res;
    }
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        } else if (o == null || this.getClass() != o.getClass()) {
            return false;
        } else {
            Region obj = (Region) o;
            return chrom.equals(getChrom()) && start.equals(obj.getStart()) && end.equals(obj.getEnd());
        }
    }
    @Override
    public String toString () {
        StringBuilder sb = new StringBuilder();
        sb.append(chrom).append("\t").append(start).append("\t").append(end);
        return sb.toString();
    }
    @Override
    public Object clone() throws CloneNotSupportedException {
        Bed obj = null;
        try {
            obj = (Bed) super.clone();
        } catch (CloneNotSupportedException ex){
            ex.printStackTrace();
        }
        return obj;
    }
    public boolean contains (String c, N pos) { // [) 左闭右开的区间
        if (chrom.equals(c)) {
            return pos.compareTo(start) >= 0 && pos.compareTo(end) < 0;
        }
        return false;
    }
    public boolean contains (String c, N s, N e){
        return contains(c,s) && contains(c,e);
    }
    public boolean contains (Bed<N,N> o) {
        return contains(o.getChrom(),o.getStart()) && contains(o.getChrom(),o.getEnd());
    }
    public boolean intersects (String c, N s, N e){
        if (chrom.equals(c)) {
            return e.compareTo(start) > 0 && s.compareTo(end) < 0;
        }
        return false;
    }
    public boolean intersects(Bed<N,N> o){
        return intersects(o.getChrom(),o.getStart(),o.getEnd());
    }

    public void mergeRegion(N s, N e){
        if (start.compareTo(s) > 0) {
            setStart(s);
        }
        if (end.compareTo(e) < 0) {
            setEnd(e);
        }
    }
    public void mergeRegion(Bed<N,N> o){
        mergeRegion(o.getStart(),o.getEnd());
    }
}



public class Region<N extends Number & Comparable,T> implements Cloneable,Comparable<Region> {
    private N start;
    private N end;
    private T value;
    public Region(N s, N e){
        start = s;
        end = e;
        if(start.compareTo(end) > 0 ){ // 序列从0开始
            throw new IllegalArgumentException("beginning of range must be less and equal than end");
        }
        value = null;
    }
    public Region(N s, N e, T data){
        start = s;
        end = e;
        if(start.compareTo(end) > 0 ){ // 序列从0开始
            throw new IllegalArgumentException("beginning of range must be less and equal than end");
        }
        value = data;
    }
    @Override
    public int compareTo(Region o) {
        int res;
        if (this == o) {
            res = 0;
        } else {
            int result = start.compareTo(o.getStart());
            if (result == 0) {
                res = end.compareTo(o.getEnd());
            } else {
                res = result;
            }
        }
        return res;
    }
    @Override
    public boolean equals(Object o) {
        boolean result;
        if (this == o) {
            result = true;
        } else if (o == null || this.getClass() != o.getClass()) {
            result = false;
        } else {
            Region obj = (Region) o;
            result = start.equals(obj.getStart()) && end.equals(obj.getEnd());
        }
        return result;
    }
    @Override
    public int hashCode () {
        int result = start.hashCode();
        result = 31*result + end.hashCode();
        result = 31*result + value.hashCode();
        return result;
    }
    @Override
    public String toString () {
        StringBuilder sb = new StringBuilder();
        sb.append(start).append("\t").append(end);
        if (value != null){
            sb.append("\t").append(value);
            // return String.format("%d\t%d",start,end);
        }
        return sb.toString();
    }
    @Override
    public Object clone() throws CloneNotSupportedException {
        Region obj = null;
        try {
            obj = (Region) super.clone();
        } catch (CloneNotSupportedException ex){
            ex.printStackTrace();
        }
        return obj;
    }
    public boolean contains (N pos) { // [) 左闭右开的区间
        return pos.compareTo(start) >= 0 && pos.compareTo(end) < 0;
    }
    public boolean contains (N s, N e){
        return contains(s) && contains(e);
    }
    public boolean contains (Region<N,N> o) {
        return contains(o.getStart()) && contains(o.getEnd());
    }
    public boolean intersects (N s, N e){
        return e.compareTo(start) > 0 && s.compareTo(end) < 0;
    }
    public boolean intersects(Region<N,N> o){
        return intersects(o.getStart(),o.getEnd());
    }
    public void setStart(N start) {
        this.start = start;
    }
    public N getStart() {
        return start;
    }
    public void setEnd(N end) {
        this.end = end;
    }
    public N getEnd() {
        return end;
    }
    public void setValue(T value) {
        this.value = value;
    }
    public T getValue() {
        return value;
    }
}
