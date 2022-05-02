import sys
from urllib import request
from scrapy.selector import Selector
import time

def deal(infile):
    sample = []
    myurl = 'https://www.ncbi.nlm.nih.gov/nuccore/'
    with open(infile,'r') as fl:
        for line in fl:
            name = line.split()
            sample.append(name[0])
    for i in sample:
        url = myurl+i
        response = request.urlopen(url)
        html = response.read()
        print(i)
        print(Selector(text=html).xpath('//title').extract())
        time.sleep(1)

if __name__ == "__main__":
    soft,infile = sys.argv
    deal(infile)
