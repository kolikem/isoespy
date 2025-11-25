import subprocess
import sys
# arg1 - r script for execution of edgeR
# arg2 - meta data
# arg3 - count data
# arg4 - output (tests)
# arg5 - output (CPM)

def de(r_script, meta_data, count_data, result_test, result_cpm):
    # Rスクリプトを実行
    subprocess.run(["Rscript", r_script, meta_data, count_data, result_test, result_cpm])

if __name__ ==  "__main__":
    de(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    print("edgeR analysis done!")
