
def read_test_file(filepath):
    mtv = []
    with open(filepath, 'r') as f:
        lines = [ l for l in f.readlines() if not l.startswith('#') ] 
        bits = int(lines[0])
        for l in lines[1:]:
            m,t,v = l.strip().split(' ')
            mtv.append((m,t,v))

    return (bits, mtv)
