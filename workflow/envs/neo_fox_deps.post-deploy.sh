#!/usr/bin/env bash
set -euo pipefail

# set all the necessary conda paths and
# ensure they exist
CONDA_BIN="${CONDA_PREFIX}/bin/"
CONDA_MAN1="${CONDA_PREFIX}/share/man/man1/"
mkdir -p $CONDA_MAN1
CONDA_INFO="${CONDA_PREFIX}/share/info/"
mkdir -p $CONDA_INFO
CONDA_LIB="${CONDA_PREFIX}/lib/"
mkdir -p ${CONDA_LIB}
CONDA_ETC="${CONDA_PREFIX}/etc/"
mkdir -p $CONDA_ETC

# install MixMHCpred, following:
# https://github.com/GfellerLab/MixMHCpred/blob/v2.1/README
MIX_MHC_PRED_VERSION="2.1"
MIX_MHC_PRED_LIB_PATH="$CONDA_PREFIX/lib/mix_mhc_pred/"
wget https://github.com/GfellerLab/MixMHCpred/archive/refs/tags/v${MIX_MHC_PRED_VERSION}.tar.gz
tar xzf v${MIX_MHC_PRED_VERSION}.tar.gz
cd MixMHCpred-${MIX_MHC_PRED_VERSION}
g++ -O3 lib/MixMHCpred.cc -o lib/MixMHCpred.x
# TODO: when updating to v2.2, change this line to:
# MMP_PLACEHOLDER="/PATH_TO_MIXMHCPRED/lib"
MMP_PLACEHOLDER="YOUR PATH TO MixMHCpred/lib FOLDER"
grep "${MMP_PLACEHOLDER}" MixMHCpred
sed -i "s%${MMP_PLACEHOLDER}%${MIX_MHC_PRED_LIB_PATH}%" MixMHCpred
mv lib $MIX_MHC_PRED_LIB_PATH
mv MixMHCpred ${CONDA_BIN}
# TODO: when updating to v2.2, change this line to:
# mv MixMHCpred_license.pdf ${CONDA_INFO}/MixMHCpred_license.pdf
mv license.pdf ${CONDA_INFO}/MixMHCpred_license.pdf
MixMHCpred -i test/test.fa -o test/out.txt -a A0101,A2501,B0801,B1801
diff <(sed '4d' test/out.txt) <(sed '4d' test/out_compare.txt)
cd ..
rm v${MIX_MHC_PRED_VERSION}.tar.gz
rm -r MixMHCpred-${MIX_MHC_PRED_VERSION}

# install MixMHC2pred, mostly following (we use the default GitHub-created
# .tar.gz file instead of the hand-crafted .zip file for future-proofing):
# https://github.com/GfellerLab/MixMHC2pred/blob/v1.2/README.md
MIX_MHC_TWO_PRED_VERSION="1.2"
MIX_MHC_TWO_PRED_LIB_PATH="${CONDA_LIB}/mix_mhc_two_pred/"
wget https://github.com/GfellerLab/MixMHC2pred/archive/refs/tags/v${MIX_MHC_TWO_PRED_VERSION}.tar.gz
tar xzf v${MIX_MHC_TWO_PRED_VERSION}.tar.gz
cd MixMHC2pred-${MIX_MHC_TWO_PRED_VERSION}
mv -t ${CONDA_BIN} MixMHC2pred MixMHC2pred_unix
mv rpep ${CONDA_ETC}
ln -s ${CONDA_ETC}/rpep ${CONDA_BIN}/rpep
mv LICENSE ${CONDA_INFO}/MixMHC2pred_unix_LICENSE
MixMHC2pred_unix -i test/testData.txt -o test/out.txt -a DRB1_11_01 DRB3_02_02 DPA1_01_03__DPB1_04_01 DQA1_05_05__DQB1_03_01
diff test/out.txt test/out_compare.txt
cd ..
rm v${MIX_MHC_TWO_PRED_VERSION}.tar.gz
rm -r MixMHC2pred-${MIX_MHC_TWO_PRED_VERSION}

# install PRIME, mostly following (minor corrections):
# https://github.com/GfellerLab/PRIME/blob/v1.0/README
PRIME_VERSION="1.0"
PRIME_LIB_PATH="${CONDA_LIB}/prime/"
wget https://github.com/GfellerLab/PRIME/archive/refs/tags/v${PRIME_VERSION}.tar.gz
tar xzf v${PRIME_VERSION}.tar.gz
cd PRIME-${PRIME_VERSION}
PRIME_PLACEHOLDER="/app/PRIME/lib"
grep "${PRIME_PLACEHOLDER}" PRIME
sed -i "s%${PRIME_PLACEHOLDER}%${PRIME_LIB_PATH}%" PRIME
mv lib $PRIME_LIB_PATH
mv PRIME ${CONDA_BIN}
mv PRIME_license.pdf ${CONDA_INFO}
PRIME -i test/test.txt -o test/out.txt -a A0201,A0101
diff <(sed '4d' test/out.txt) <(sed '4d' test/out_compare.txt)
cd ..
rm v${PRIME_VERSION}.tar.gz
rm -r PRIME-${PRIME_VERSION}

# This is the non-portable version of the 1st line
# in both netMHCpan and netMHCIIpan, assuming a
# root install of tcsh. 
TCSH_ROOT="#! /bin/tcsh -f"
# For the scripts to work with any tcsh in the
# $PATH, we need to change this to:
TCSH_PATH="#!/usr/bin/env tcsh"

# install netMHCpan version 4.1
# requires tcsh to have been installed via conda dependencies
NET_MHC_PAN_4_1_LIB="${CONDA_LIB}/netMHCpan_4_1/"
mkdir -p ${NET_MHC_PAN_4_1_LIB}
NET_MHC_PAN_4_1_ETC="${CONDA_ETC}/netMHCpan_4_1/"
mkdir -p ${NET_MHC_PAN_4_1_ETC}
tar xzf ${NET_MHC_PAN_4_1_TARBALL}
cd netMHCpan-4.1
wget https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz
tar xzf data.tar.gz
rm data.tar.gz
grep "${TCSH_ROOT}" netMHCpan
sed -i "s%${TCSH_ROOT}%${TCSH_PATH}%" netMHCpan
grep -P "setenv\s+NMHOME" netMHCpan
sed -r -i "s%^setenv\s+NMHOME+.*$%setenv NMHOME ${NET_MHC_PAN_4_1_LIB}%" netMHCpan
mv -t ${CONDA_BIN} netMHCpan
mv -t ${CONDA_MAN1} netMHCpan.1
mv -t ${CONDA_INFO} netMHCpan-4.1.readme
mv -t ${NET_MHC_PAN_4_1_LIB} Linux_x86_64
mv -t ${NET_MHC_PAN_4_1_ETC} data
ln -s ${NET_MHC_PAN_4_1_ETC}/data ${NET_MHC_PAN_4_1_LIB}/data
cd test
netMHCpan -p test.pep -BA -xls -a HLA-A01:01,HLA-A02:01 -xlsfile my_NetMHCpan_out.xls
diff NetMHCpan_out.xls my_NetMHCpan_out.xls
cd ../..
rm -r netMHCpan-4.1

# install netMHCIIpan version 4.0
# requires tcsh to have been installed via conda dependencies
NET_MHC_TWO_PAN_4_0_LIB="${CONDA_LIB}/netMHCIIpan_4_0/"
mkdir -p ${NET_MHC_TWO_PAN_4_0_LIB}
NET_MHC_TWO_PAN_4_0_ETC="${CONDA_ETC}/netMHCIIpan_4_0/"
mkdir -p ${NET_MHC_TWO_PAN_4_0_ETC}
tar xzf ${NET_MHC_TWO_PAN_4_0_TARBALL}
cd netMHCIIpan-4.0
wget https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/data.tar.gz 
tar xzf data.tar.gz
rm data.tar.gz
grep "${TCSH_ROOT}" netMHCIIpan
sed -i "s%${TCSH_ROOT}%${TCSH_PATH}%" netMHCIIpan
grep -P "setenv\s+NMHOME" netMHCIIpan
sed -r -i "s%^setenv\s+NMHOME+.*$%setenv NMHOME ${NET_MHC_TWO_PAN_4_0_LIB}%" netMHCIIpan
mv -t ${CONDA_BIN} netMHCIIpan
mv -t ${CONDA_MAN1} netMHCIIpan.1
mv -t ${CONDA_INFO} netMHCIIpan-4.0.readme
mv -t ${NET_MHC_TWO_PAN_4_0_LIB} Linux_x86_64 NetMHCIIpan-4.0.pl
mv -t ${NET_MHC_TWO_PAN_4_0_ETC} data
ln -s ${NET_MHC_TWO_PAN_4_0_ETC}/data ${NET_MHC_TWO_PAN_4_0_LIB}/data
cd test
netMHCIIpan -f example.fsa -a DRB1_0101 > example.fsa.myout
diff example.fsa.out example.fsa.myout
netMHCIIpan -f example.pep -inptype 1 -a DRB1_0101 > example.pep.myout
diff example.pep.out example.pep.myout
netMHCIIpan -f example.fsa -a H-2-IAb -s -u > example.fsa.sorted.myout
diff example.fsa.sorted.out example.fsa.sorted.myout
netMHCIIpan -f example.fsa -hlaseq DRB10101.fsa > example.fsa_hlaseq.myout
diff example.fsa_hlaseq.out example.fsa_hlaseq.myout
netMHCIIpan -f example.fsa -hlaseqA alpha.dat -hlaseq beta.dat > example.fsa_hlaseq_A+B.myout
diff example.fsa_hlaseq_A+B.out example.fsa_hlaseq_A+B.myout
cd ../..
rm -r netMHCIIpan-4.0
