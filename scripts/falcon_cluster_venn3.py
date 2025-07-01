import sys
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

map_set_file = {}
lines = open(sys.argv[1], 'r').readlines()
for l in lines:
  es = l.strip().split()
  map_set_file[ es[0] ] = es[1]
#print( map_set_file )

lines = [ l for l in open(sys.argv[2], 'r').readlines() if not l.startswith('#') ]
t = {}
for i, e in enumerate(lines[0].strip().split(',')):
  t[e] = i
#print(t)

map_cluster_files = defaultdict(list)
for l in lines[1:]:
  es = l.strip().split(',')
  ident = es[t['identifier']]
  fname = ident.split(':')[2]
  cluster = int(es[t['cluster']])
  if cluster > -1 and fname in map_set_file.keys():
    map_cluster_files[cluster].append(map_set_file[fname])

#print(map_cluster_files)

# 定义各区域的数值
# 格式: (A独有, B独有, AB交集, C独有, AC交集, BC交集, ABC交集)
counts = {
    '100': 0,  # 仅A
    '010': 0,  # 仅B
    '110': 0,   # AB交集
    '001': 0,  # 仅C
    '101': 0,   # AC交集
    '011': 0,   # BC交集
    '111': 0    # ABC交集
}

# 获取三个文件的标识（假设map_set_file.values()返回的是[A, B, C]）
file_ids = list(map_set_file.values())
if len(file_ids) != 3:
    raise ValueError("需要正好三个文件标识")

A, B, C = file_ids[0], file_ids[1], file_ids[2]

# 统计
for cluster, files in map_cluster_files.items():
    files = set(files)
    if len(files) == 1:
        if A in files:
            counts['100'] += 1
        elif B in files:
            counts['010'] += 1
        elif C in files:
            counts['001'] += 1
    elif len(files) == 2:
        if A in files and B in files and C not in files:
            counts['110'] += 1
        elif A in files and C in files and B not in files:
            counts['101'] += 1
        elif B in files and C in files and A not in files:
            counts['011'] += 1
    elif len(files) == 3:
        counts['111'] += 1

# 将字典转换为venn3需要的元组格式
venn_values = (
    counts['100'], counts['010'], counts['110'],
    counts['001'], counts['101'], counts['011'],
    counts['111']
)

# 创建Venn图
venn = venn3(subsets=venn_values, set_labels=map_set_file.values())

# 设置所有文本的字体大小（包括区域数字和集合标签）
for text in venn.set_labels:
    text.set_fontsize(14)  # 增大集合标签字体大小（A/B/C）
for text in venn.subset_labels:
    text.set_fontsize(12)  # 增大区域数字字体大小

# 添加圆形轮廓
venn3_circles(subsets=venn_values)

# 显示图形
plt.title("Venn Diagram of Three Sets")
plt.show()
