l = list()
l2 = list()
seq = str()
with open('./test.file', 'r') as f:
    for line in f.readlines():
        seq += line
    print([seq])
    seq = seq.replace('\n', '')
    print([seq])

