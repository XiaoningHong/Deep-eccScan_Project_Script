import os,sys,re,gzip,csv
import pandas as pd
from argparse import ArgumentParser
import subprocess
import pandas as pd

from typing import Dict, Tuple


def analyze_coords_file(coords_file: str, output_file: str, similarity_threshold: float = 90.0) -> Dict:
    """
    分析nucmer比对坐标文件
    
    参数:
    coords_file: 坐标文件路径
    similarity_threshold: 相似度阈值（百分比）
    
    返回:
    包含分析结果的字典
    """
    if not os.path.exists(coords_file):
        logger.error(f"坐标文件不存在: {coords_file}")
        return {}
    bed1_count = set(); bed2_count = set();
    BP_diff = 0;count = 0;
    try:
        with open(coords_file, 'r') as f, open(output_file, 'w') as f_output:
            lines = f.readlines()
            # [S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [TAGS]
            for line in lines:
                parts = re.split("\t", line.strip())
                similarity = float(parts[6])  
                if similarity >= similarity_threshold and float(parts[9]) >= similarity_threshold and float(parts[10]) >= similarity_threshold and (int(parts[7])+int(parts[8])-int(parts[4])-int(parts[5])<250):
                    f_output.write(line)
                    count += 1
                    BP_diff += (int(parts[7]) - int(parts[4])) + (int(parts[8]) - int(parts[5])) #(LENR - LEN1) + (LENQ - LEN2)
                    bed1_count.add(parts[-2]); bed2_count.add(parts[-1]);
            mean_BP_diff = round(BP_diff / count, 2)
            analysis_results = {
            'file1_matches': bed1_count,
            'file2_matches': bed2_count,
            'Base_pair_Difference': mean_BP_diff
            }
        
        return analysis_results
        
    except Exception as e:
        analysis_results = {
            'file1_matches': bed1_count,
            'file2_matches': bed2_count,
            'Base_pair_Difference': 'NULL'
            }
        return analysis_results

def read_bed(bed_file: str):
    bed_dict = {}
    with open(bed_file, 'r') as f_:
        for line in f_:
            ls  = re.split("\t", line.strip())
            _id = f"{ls[0]}:{ls[1]}-{ls[2]}"
            bed_dict[_id] = ls

    return bed_dict

def evaluate_model(TP, TP2, FP, TN, FN):
    """模型评估"""
    total = TP + FP + TN + FN
    
    # 安全除法函数
    def safe_div(a, b):
        return a / b if b != 0 else 0.0
    
    accuracy = safe_div(TP + TN, total)
    precision = safe_div(TP, TP + FP)
    recall = safe_div(TP, TP + FN)
    FPR = safe_div(FP, FP + TN)
    specificity = safe_div(TN, TN + FP)
    FNR = safe_div(FN, TP +FN)
    f1_score = safe_div(2 * precision * recall, precision + recall)
    
    print("评估结果:")
    print(f"Accuracy: {accuracy:.4f} ({accuracy*100:.2f}%)")
    print(f"Precision: {precision:.4f} ({precision*100:.2f}%)")
    print(f"Recall: {recall:.4f} ({recall*100:.2f}%)")
    print(f"False Positive Rate: {FPR:.4f} ({FPR*100:.2f}%)")
    print(f"Specificity: {specificity:.4f} ({specificity*100:.2f}%)")
    print(f"False Negative Rate: {FNR:.4f} ({FNR*100:.2f}%)")
    print(f"F1 score: {f1_score:.4f}")
    
    return {
        'TP1': TP, 'TP2': TP2, 'FP':FP, 'TN':TN, 'FN':FN,
        'accuracy': accuracy,
        'precision': precision, 'false_positive_rate': FPR,
        'recall': recall,
        'specificity': specificity, 'false_negative_rate' : FNR,
        'f1_score': f1_score
    }

def save_bed(regions, filename, bed_dict=None):
    """BED保存函数，如果region在bed_dict中，则输出字典中的value"""
    
    data = []
    for region in regions:
        chrom, positions = region.split(':')
        start, end = positions.split('-')
        
        # 准备行数据
        if bed_dict is not None and region in bed_dict:
            data.append((chrom, int(start), int(end), "\t".join(bed_dict[region][3:]), region))
        else:
            data.append((chrom, int(start), int(end), None, region))
    
    # 按chrom和start排序
    data.sort(key=lambda x: (x[0], x[1]))
    
    with open(filename, 'w') as f:
        for chrom, start, end, info, region in data:
            if info is not None:
                line = f"{chrom}\t{start}\t{end}\t{info}"
            else:
                line = f"{chrom}\t{start}\t{end}"
            f.write(line + "\n")

def get_FP_Number(regions, match_bed):
    """
    Count regions with no overlap in match_bed
    """
    count = {}
    
    for region in regions:
        chrom, pos = region.split(':')
        r_start, r_end = map(int, pos.split('-'))
        
        # Check if overlapped by any bed region
        for b_region in match_bed:
            b_chrom, b_pos = b_region.split(':')
            b_start, b_end = map(int, b_pos.split('-'))
            b_start, b_end = int(b_start), int(b_end)
            
            if (chrom == b_chrom and 
                r_end > b_start and r_start < b_end):
                break  # Found overlap
        else:
            # No overlap found for this region
            count[region]= 1
    
    return count
        
    
def main(args):
    bed1 = args.bed1; bed2 = args.bed2; bed3 = args.bed3;
    tolerance = args.threshold; 
    threads = args.threads; 
    output_dir = args.output;
    
    #Bedtools
    subprocess.run(f"/opt/anaconda3/bin/bedtools getfasta -fi /home/DataBase/Refrence/hg38/hg38.fa -bed {bed1} -fo {output_dir}/Positive.fasta", shell=True)
    subprocess.run(f"/opt/anaconda3/bin/bedtools getfasta -fi /home/DataBase/Refrence/hg38/hg38.fa -bed {bed2} -fo {output_dir}/Negative.fasta", shell=True)
    subprocess.run(f"/opt/anaconda3/bin/bedtools getfasta -fi /home/DataBase/Refrence/hg38/hg38.fa -bed {bed3} -fo {output_dir}/Detection.fasta", shell=True)
    #Mummer Positive
    subprocess.run(f"/home/DataBase/miniforge3/envs/mummer4/bin/nucmer --maxmatch -l 100 -c 100  -t {threads} -p {output_dir}/mummer.pos {output_dir}/Positive.fasta {output_dir}/Detection.fasta", shell=True)  
    subprocess.run(f"/home/DataBase/miniforge3/envs/mummer4/bin/show-coords -rclTH {output_dir}/mummer.pos.delta > {output_dir}/mummer.pos.coords", shell=True)
    #Mummer Negative
    subprocess.run(f"/home/DataBase/miniforge3/envs/mummer4/bin/nucmer --mum -l 100 -c 100  -t {threads} -p {output_dir}/mummer.neg {output_dir}/Negative.fasta {output_dir}/Detection.fasta", shell=True)  
    subprocess.run(f"/home/DataBase/miniforge3/envs/mummer4/bin/show-coords -rclTH {output_dir}/mummer.neg.delta > {output_dir}/mummer.neg.coords", shell=True)
    
    
    print("=" * 50)
    print("Positive BED文件比较结果")
    print("=" * 50)
    result_pos = analyze_coords_file(f"{output_dir}/mummer.pos.coords", f"{output_dir}/mummer.pos.matched.coords")
    result_neg = analyze_coords_file(f"{output_dir}/mummer.neg.coords", f"{output_dir}/mummer.neg.matched.coords")

    pos_set = set(read_bed(bed1))
    neg_set = set(read_bed(bed2))
    result_set = set(read_bed(bed3))

    
    print(f"Positive 总区域数: {len(pos_set)}")
    print(f"检测结果总区域数: {len(result_set)}")
    #Positive
    TP1 = len(result_pos['file1_matches']); TP1_ratio = TP1/len(pos_set)
    pos_unmatched = pos_set - result_pos['file1_matches']
    print (f"TP1: {TP1}\nMummer:检测结果中真阳性的比例: {TP1_ratio:.2%}")
    FN = len(pos_unmatched) 
    FN_ratio = FN / len(pos_set) if FN > 0 else 0.0
    print(f"检测结果中不匹配 Positive 区域数: {FN}")
    print (f"FN: {FN}\nMummer:检测结果中假阴性比例: {FN_ratio:.2%}")
    save_bed(pos_unmatched, f"{output_dir}/Positive_unmatched.bed")
    #Software
    TP2 = len(result_pos['file2_matches']); TP2_ratio = TP2/len(result_set)
    save_bed(result_pos['file2_matches'], f"{output_dir}/Software_matched.bed", read_bed(bed3))
    result_unmatched = result_set - result_pos['file2_matches']
    
    FP_set = set(get_FP_Number(result_unmatched, result_pos['file2_matches']))
    FP = len(FP_set)
    FP = len(result_unmatched)
    print (f"TP2: {TP2}\nMummer:检测结果中 Positive 匹配比例: {TP2_ratio:.2%}")
    print(f"检测结果不匹配 Positive 区域数: {FP}\nFP: {FP}")
    save_bed(FP_set, f"{output_dir}/Software_unmatched.bed", read_bed(bed3))
    #Negative
    TN = len(neg_set - (result_neg['file2_matches']- result_set));
    TN_ratio = TN/len(neg_set) if TN > 0 else 0.0
    print(f"检测结果不匹配 Nagetive 区域数: {TN}\t比例: {TN_ratio:.2%}")
    print("=" * 50)
    print(f"TP: {TP2}\nFP: {FP}\nTN: {TN}\nFN: {FN}\n")
    metrics = evaluate_model(TP1, TP2, FP, TN, FN)
    #
    metrics['Base_pair_Difference'] = result_pos['Base_pair_Difference']
    print(f"Base pair Difference: {result_pos['Base_pair_Difference']}")
    DupRatio = round((TP2/TP1), 2)
    print(f"Redundancy: {DupRatio}")
    metrics['Redundancy'] = DupRatio
    df = pd.DataFrame(list(metrics.items()), columns=['metric', 'value'])
    df.to_csv(f"{output_dir}/Metrics.tsv", sep='\t', index=False)
    print("=" * 50)

if __name__ == "__main__":
    parser = ArgumentParser(description='Example for python script')
    parser.add_argument("-1", "--bed1", action="store", dest="bed1",
                        help="Input Bed 1 of Positive", required=True)
    parser.add_argument("-2", "--bed2", action="store", dest="bed2",
                        help="Input Bed 2 of Negative", required=True)
    parser.add_argument("-3", "--bed3", action="store", dest="bed3",
                        help="Input Bed 3 of Software", required=True)
    parser.add_argument("-o", "--output", action="store", dest="output",
                        help="output", required=True)
    parser.add_argument('-c', '--threshold', type=float, default=50, 
                       help='Tolerance (default:50)')
    parser.add_argument("-t", "--threads", type=int, dest="threads", default=64,
                        help=f"Number of processes to use (default: 64)")

    args = parser.parse_args()

    main(args)
