import urllib.request
from bs4 import BeautifulSoup
import urllib
import re
from Bio import SeqIO
from Bio import Seq


def main():

    url = "file:///home/wanghm/Desktop/BlockLogo.html"

    html = urllib.request.urlopen(url)

    bsObj = BeautifulSoup(html,features='html.parser')

    nuc_numer = bsObj.find_all(name="td", text=re.compile(r".*A.*"))

    with open("/home/wanghm/Desktop/primary_motif.fasta", "w+") as f:
        for idx,i in enumerate(nuc_numer):
            if idx>1 and idx<=138:
                ss=SeqIO.SeqRecord(Seq.Seq(i.get_text()),id="%s"%idx, description="")
                SeqIO.write(ss, f, format="fasta")

if __name__ == '__main__':
    main()
