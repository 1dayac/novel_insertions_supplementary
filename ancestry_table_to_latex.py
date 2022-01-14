with open("latex.txt", 'w') as latex:
    with open("ancestry_specific_insertions.txt", 'r') as table:
        for r in table.readlines():
            split = r.strip().split()
            race = split[0]
            pos = split[1].split("_")[0] + " " + split[6]
            num_1 = split[2]
            num_2 = split[3]
            freq = str(round(float(split[5]), 4))
            line = " & ".join([race, pos, num_1, num_2, freq])
            latex.write(line + "\\\\\n" )