f = open("mult_result.txt","r")
test = []
for line in f:
    l = line.strip()
    test.append(l)

for r in test:
    print r[0]
    print r[1]
    print r[2]
