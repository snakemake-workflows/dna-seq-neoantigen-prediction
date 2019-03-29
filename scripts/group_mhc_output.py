import sys
import os

def group_mhc_output(mhcout):
    al = next(mhcout).split('\t')
    header = next(mhcout).rstrip().split('\t')
    allele = ""
    for i in range(0,len(header)):
        if i < len(al)-1:
            if al[i] != '':
                allele = "_" + al[i]
            header[i] += allele

    print('\t'.join(header))

    for line in mhcout:
        print(line.rstrip())

def main():
    if os.path.getsize(sys.argv[1]) > 0:
        mhcout = open(sys.argv[1],'r')
        group_mhc_output(mhcout)
    else:
        print("")

if __name__ == '__main__':
    sys.exit(main())
