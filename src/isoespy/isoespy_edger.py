import argparse
import subprocess
try:
    from run_edgeR_n_get_results import de
except ModuleNotFoundError:
    from isoespy.run_edgeR_n_get_results import de

def main(args=None):
    # parser
    parser = argparse.ArgumentParser(description='isoespy_edger')
    parser.add_argument('-r', '--r_script', required=True, type=str, default=None, help='R script for edgeR')
    parser.add_argument('-m', '--meta_data', required=True, type=str, default=None, help='metadata')
    parser.add_argument('-c', '--count_data', required=True, type=str, default=None, help='count expression data')
    parser.add_argument('-o1', '--output_de', required=True, type=str, default=None, help='result of edgeR DE analysis')
    parser.add_argument('-o2', '--output_cpm', required=True, type=str, default=None, help='result of edgeR CPM expression data')

    args = parser.parse_args()
    r_script = args.r_script
    meta_data = args.meta_data
    count_data = args.count_data
    det_data = args.output_de
    cpm_data = args.output_cpm

    # DE analysis
    de(r_script, meta_data, count_data, det_data, cpm_data)

    # adjustment
    det_data_reorder = det_data + ".reorder"
    subprocess.run([
        "awk", "-F,",
        '{print $6","$1","$5","$2","$3","$4}',
        det_data
        ], stdout=open(det_data_reorder, "w"))
    subprocess.run(["awk", "NR==1{$0=\"#\"$0}1", det_data_reorder], stdout=open("temp.txt", "w"))
    subprocess.run(["mv", "temp.txt", det_data_reorder])


if __name__ == "__main__":
    main()

