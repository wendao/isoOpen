import sys
from collections import defaultdict

lines = open(sys.argv[1], 'r').readlines()
t = {}
groups = defaultdict(lambda: defaultdict(list))
for l in lines:
  l = l.strip()
  if len(l) == 0: continue
  if l[0] == "#": continue

  if not t:
    for i, k in enumerate(l.split(',')):
      t[k] = i
    continue

  es = l.split(',')
  cluster = int(es[t['cluster']])
  if cluster == -1: continue

  identifier = es[t['identifier']].replace("mzspec:USI000000:20250606_PKA_", "")
  #print( identifier, cluster )
  fn, _, scan = identifier.split(':')
  groups[cluster][fn].append(int(scan))

for clust in groups.keys():
  if len(groups[clust].keys())==5:
    #print(groups[clust])
    out_str = "cluster:"
    for k in groups[clust].keys():
      out_str += " "+k
      for s in groups[clust][k]:
        out_str += ":"+str(s)
    print(out_str)
