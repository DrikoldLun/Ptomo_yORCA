lines = []
with open('Z_mag_6.0_10.0.dat','r') as f:
    for line in f.readlines():
        if len(line.split()) < 16:
            continue
        lines.append(line)
    f.close()
with open('Z_mag_6.0_10.0.dat','w+') as f:
    f.writelines(lines)
    f.close()

