import numpy as np

gap.load_package("grape")
gap.Read('"ziv-av.g"')

__author__ = 'tony'

input = sys.argv[1]
output = sys.argv[2]
output2 = sys.argv[3]

s = '"' + input + '"'
gap.Read(s)

fw = open(output, 'w')
fw2 = open(output2, 'w')

size = gap('Length(sol);')
rang = list(range(1, size + 1))
comp = []
rozklad = []

for x in range(1, size+1):

    if x in rang:

        if len(rang) != 0:
            print x, len(rang)
        else:
            break

        newrang = set()


        for y in rang:

            pom = gap.IsIsomorphicCgr_2(gap('List(sol['+str(x)+'],ShallowCopy)'), gap('List(sol['+str(y)+'],ShallowCopy)'))

            if pom:
                newrang.add(y)

                comp.append(gap('sol['+str(y)+']'))

        rang = [x for x in rang if x not in newrang]

        if comp != []:
            rozklad.append(comp[0])
            comp = []

print "-----------------------"
print "number:", len(rozklad)

fw2.write(str(len(rozklad)))

st = ""
st += ("sol:=[")
for i in rozklad:

    st += str(i) + ","

st = st[:-1]
st += "];"
fw.write(st)

fw.close()
fw2.close()
