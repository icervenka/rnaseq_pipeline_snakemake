def print_header(s):
    header_len = len(s)
    print("╔"+78*"═"+"╗")
    print("║"+"  "+s+"  "+(74-header_len)*" "+"║")
    print("╚"+78*"═"+"╝")


def all_paired(Metadata):
    r1 = sum(Metadata["read"] == "R1")
    r2 = sum(Metadata["read"] == "R2")
    if r1 == r2:
        return True
    else:
        return False
    

def all_single(Metadata):
    r2 = sum(Metadata["read"] == "R2")
    if r2 == 0:
        return True
    else:
        return False

        
def is_set_subsample(s):
    if s == "" or s == 0 or s == None:
        return False
    elif s in ['False', 'false', 'FALSE']:
        return False
    else:
        return True


def is_set_trimmer(s):
    if s in ['fastp', 'trimmomatic', 'cutadapt']:
        return True
    else:
        return False