import sys

exp = []
obs = []
for line in sys.stdin:
    if ">" in line:
        if tuple(line.replace(">", "").split()[0].split("_")) not in exp:
            exp.append(tuple(line.replace(">", "").split()[0].split("_")))
    else:
        if tuple(line.split()[:2]) not in obs:
            obs.append(tuple(line.split()[:2]))

print "Se="+str(len(list(set(exp).intersection(obs)))/float(len(exp))),
print "Sp="+str(len(list(set(exp).intersection(obs)))/float(len(obs)))
