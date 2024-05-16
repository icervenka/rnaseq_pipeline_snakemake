def print_header(s):
    header_len = len(s)
    print("╔"+78*"═"+"╗")
    print("║"+"  "+s+"  "+(74-header_len)*" "+"║")
    print("╚"+78*"═"+"╝")

def has_extra_config(conf, conf_extra):
    conf = conf.strip()
    if len(conf) > 0:
        return conf_extra[conf]
    else:
        return " "