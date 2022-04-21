#!/public/agis/huangsanwen_group/wangxin/bin/python3

import sys
n=int(sys.argv[1])
spe=sys.argv[2]

for i in range(1,n+1):
    print('/vol1/agis/huangsanwen_group/lihongbo/software/GlimmerHMM/bin/glimmerhmm_linux_x86_64 ./seq/'+str(i)+'.fa -d cucumber_ok  -g -o ./glimmerhmm/'+str(i)+'.out  # 01_train  cucumber_ok directory')
