import json
import os,sys,time

def get_matrix_data(overlap_json):
    starttime = time.time()
    with open(overlap_json, 'r') as fm:
        matrix_data = json.load(fm)
    endtime = time.time()
    sys.stdout.write("Loading TR matrix file: {} seconds \n ".format(endtime-starttime))
    return matrix_data


overlap_json='/project/zanglab_project/hz9fq/annotations/bart_library/mm10_library/bart2_mm10_TF_overlap.json'
overlap_dict = get_matrix_data(overlap_json)

lstring="".join(overlap_dict.values())
tf_occurance_str=lstring.strip().split(' ')

tf_occurance=[int(x) for x in tf_occurance_str]

del tf_occurance_str

tf_len = max(tf_occurance)
print(tf_len)

#sort the occurance list
tf_occurance.sort()

tf_binding_count={}

for i in range(1,tf_len):
    count=0
    while True:
        if tf_occurance[count]==i:
            count+=1
            continue
        else:
            print(i,count)
            tf_binding_count[i]=count
            del tf_occurance[0:count]
            break

print(tf_len,len(tf_occurance))
tf_binding_count[tf_len]=len(tf_occurance)
 
with open("bart2_mm10_binding_count.json", "w") as outfile:
    json.dump(tf_binding_count, outfile)

