import sys
from tqdm import tqdm
import collections
from pyteomics import mzml

def find_frequent_ions(
    mzml_file_path,
    target_scan_numbers,
    frequency_threshold=0.5,
    mz_tolerance=0.01
):
    """
    在指定的MS2谱图列表中查找并统计频繁出现的离子。

    参数:
    mzml_file_path (str): mzML文件的路径。
    target_scan_numbers (list): 包含目标MS2谱图Scan Number的列表。
    frequency_threshold (float): 离子出现的频率阈值 (0.0到1.0之间)。
    mz_tolerance (float): m/z值的容差（单位：Da），用于离子分箱。
    """

    # 将列表转换为集合，以提高查找效率
    target_scans_set = set(target_scan_numbers)
    total_spectra_to_analyze = len(target_scans_set)

    if total_spectra_to_analyze == 0:
        print("错误：提供的Scan Number列表为空。")
        return

    print(f"准备分析 {total_spectra_to_analyze} 张指定的MS2谱图...")
    print(f"文件路径: {mzml_file_path}")
    print(f"频率阈值: > {frequency_threshold:.0%}")
    print(f"m/z 容差: +/- {mz_tolerance} Da\n")

    # 用于统计每个分箱后的离子在多少张谱图中出现过
    # collections.Counter 是一个专门用于计数的字典
    ion_occurrence_counter = collections.Counter()

    #多少张MS2谱图
    total_ms2 = max(target_scan_numbers)
    
    # 使用 pyteomics 读取 mzML 文件
    processed = 0
    with mzml.read(mzml_file_path) as reader:
        for spectrum in tqdm(reader, total=total_ms2, desc="Processing spectra"):
            # 检查是否是MS2谱图
            if spectrum.get('ms level') != 2:
                continue

            # 从谱图的ID中提取scan number
            # 常见的ID格式为 'controllerType=0 controllerNumber=1 scan=12345'
            try:
                # 这种解析方式更稳健
                scan_num_str = spectrum['id'].split('scan=')[-1]
                current_scan_num = int(scan_num_str)
                if current_scan_num > total_ms2: break
            except (KeyError, IndexError, ValueError):
                # 如果ID格式不同，跳过此谱图并给出提示
                # print(f"警告：无法从ID '{spectrum.get('id', '')}' 中解析Scan Number。")
                continue

            # 检查当前谱图是否是我们关心的目标
            if current_scan_num in target_scans_set:
                
                # 用于存放当前谱图出现过的所有离子（分箱后），使用集合避免重复计数
                ions_in_this_spectrum = set()
                
                mzs = spectrum['m/z array']
                intensities = spectrum['intensity array']

                for mz in mzs:
                    # **核心步骤：离子分箱**
                    # 通过将m/z值除以容差、取整，再乘回容差，实现分箱
                    # 例如，容差为0.01，则147.112和147.113都会被归入147.11这个箱子
                    binned_mz = round(mz / mz_tolerance) * mz_tolerance
                    ions_in_this_spectrum.add(binned_mz)
                
                # 更新全局计数器
                # Counter.update()会为集合中的每个元素计数加一
                ion_occurrence_counter.update(ions_in_this_spectrum)

    print("分析完成，正在整理统计结果...\n")

    # 整理并筛选结果
    frequent_ions = []
    for mz_bin, count in ion_occurrence_counter.items():
        frequency = count / total_spectra_to_analyze
        if frequency >= frequency_threshold:
            frequent_ions.append({
                'mz_bin': mz_bin,
                'count': count,
                'frequency': frequency
            })

    # 按频率降序排序
    frequent_ions.sort(key=lambda x: x['frequency'], reverse=True)

    # 打印结果
    print("------ 频繁出现的离子统计结果 ------")
    if not frequent_ions:
        print(f"在指定的 {total_spectra_to_analyze} 张谱图中，没有发现出现频率超过 {frequency_threshold:.0%} 的离子。")
        print("您可以尝试降低'frequency_threshold'或增大'mz_tolerance'。")
        return

    print(f"{'m/z (Binned)':<15} {'出现次数':<12} {'总谱图数':<12} {'出现频率':<10}")
    print("-" * 55)
    for ion in frequent_ions:
        # 使用 f-string 格式化输出，使表格对齐
        print(
            f"{ion['mz_bin']:<15.4f} "
            f"{ion['count']:<12} "
            f"{total_spectra_to_analyze:<12} "
            f"{ion['frequency']:.2%}"
        )

# --- 使用示例 ---
if __name__ == '__main__':
    # 1. 设置您的文件路径
    MZML_FILE = sys.argv[1]

    # 2. 提供您关心的Scan Number列表
    SCAN_NUMBERS = []
    lines = open(sys.argv[2], 'r').readlines()
    for l in lines:
        es = l.strip().split()
        SCAN_NUMBERS.append(int(es[0]))

    # 3. 设置参数
    FREQ_THRESHOLD = 0.5  # 频率阈值设为50%
    MZ_TOLERANCE = 0.01   # 质量容差设为0.01 Da

    # 运行分析函数
    find_frequent_ions(MZML_FILE, SCAN_NUMBERS, FREQ_THRESHOLD, MZ_TOLERANCE)
