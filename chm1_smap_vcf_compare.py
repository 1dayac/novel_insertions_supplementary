from common import Insertion

count = 0
with open("CHM1_bionano.smap", "r") as bionano:
    for r in bionano.readlines():
        if r.startswith("#") or r.find("deletion") != -1 or r.find("inversion") != -1 or r.find("trans") != -1:
            continue

        splited = r.split('\t')
        significance = float(splited[8])
        if significance == -1.0:
            continue

        chromID = str(splited[2])
        qstart = int(float(splited[4]))
        qend = int(float(splited[5]))
        rstart = int(float(splited[6]))
        rend = int(float(splited[7]))
        count += 1
        if abs(rstart - rend) < 5000 and abs(qend - qstart) >= 2000:
            print(str(abs(qstart - qend)) + " " + str(abs(rstart - rend)))

    print(count)