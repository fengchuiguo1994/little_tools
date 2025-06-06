#!/usr/bin/env python3
# coding: utf-8
"""
@file: spatial_cluster.py
@description: 
@author: Ping Qiu
@email: qiuping1@genomics.cn
@last modified by: JunHao Xu

change log:
    2022/01/05  create file.
    2022/03/18  change clustering method to leiden, change file name from cell_cluster.py to spatial_cluster.py.
    2022/08/11  increase scale function.
    2022/10/20  increase error code.
    2022/11/23  set pca function param use_highly_genes to False.
"""

###### Version and Date
PROG_VERSION = '1.4.0'
PROG_DATE = '2022-11-23'

import stereo as st
import argparse
import time
import os, sys
from scipy.sparse import issparse
from optparse import OptionParser

import faulthandler
faulthandler.enable()

def spatial_cluster(gef_file, bin_size):
    """
    run the cluster using stereopy.

    :param gef_file: gef file path.
    :param bin_size: bin size.
    :return: StereoExpData.
    """
    
    data = st.io.read_gef(gef_file, bin_size=bin_size)
    data.tl.cal_qc()
    data.tl.filter_cells(min_gene=1)
    data.tl.raw_checkpoint()
    # if issparse(data.exp_matrix):
    #     data.exp_matrix = data.exp_matrix.toarray()
    data.tl.normalize_total()
    data.tl.log1p()
    if len(data.gene_names) < 3000:
        print_err("Gene number less than 3000, please check your gef file","SAW-A70005")
    data.tl.highly_variable_genes(min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=3000, res_key='highly_variable_genes')
    data.tl.scale()
    data.tl.pca(use_highly_genes=False, hvg_res_key='highly_variable_genes', n_pcs=20, res_key='pca')
    data.tl.neighbors(pca_res_key='pca', n_pcs=30, res_key='neighbors')
    data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap')
    data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden')
    return data


def stereo2anndata(data, out_file=None):
    """
    transform the StereoExpData object into Anndata object.

    :param data: StereoExpData
    :param out_file: h5ad output path of Anndata object
    :return:
    """
    return st.io.stereo_to_anndata(data, output=out_file)


def print_err(case, code):
    """
    print error code
    """
    err_code={
        "SAW-A70001":"{} is missing.",
        "SAW-A70002":"cannot access {}: No such file or directory.",
        "SAW-A70003":"{} file format error.",
        "SAW-A70004":"info loss:{}.",
        "SAW-A70005":"{}."
    }
    nowtime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    with open("errcode.log","a") as err_info:
        err_info.write("[{}] {}: {}\n".format(nowtime, code, err_code[code].format(case)))
    sys.stderr.write("[{}] {}: {}\n".format(nowtime, code, err_code[code].format(case)))
    sys.exit(1)


def main():
    """
    %prog [options]
    generate cluster analysis
    Usage: python ./spatial_cluster.py -i ./SN.gef -o ./cluster.h5ad .
    """

    parser = OptionParser(main.__doc__)

    parser.add_option('-s', '--bin_size', action='store', type=int, default=50, help='The bin size or max bin szie that to combine the dnbs. default=50')
    parser.add_option('-i', '--gef_file', action='store', help='gef file path.')
    parser.add_option('-o', '--out_file', action='store', help='outfile(h5ad)')

    opts, args = parser.parse_args()

    if(opts.bin_size == None):
        print_err("-s or --bin_size","SAW-A70001")

    if(opts.gef_file == None):
        print_err("-i or --gef_file","SAW-A70001")
    elif not os.path.exists(opts.gef_file):
        print_err(os.path.abspath(opts.gef_file),"SAW-A70002")

    if(opts.out_file == None):
        print_err("-o or --out_file","SAW-A70001")

    if opts.bin_size not in [1, 10, 20, 50, 100, 200, 500]:
        print_err("The bin size is out of range, please check. The range of gef binsize is [1,10,20,50,100,200,500]","SAW-A70005")

    stereo_data = spatial_cluster(opts.gef_file, opts.bin_size)
    stereo2anndata(stereo_data, opts.out_file)


if __name__ == '__main__':
    main()
