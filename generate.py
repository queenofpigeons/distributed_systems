#generate a random nxn matrix

import random

f = open('input.txt', 'w')
n = int(input())
f.write(str(n))
f.write(" ")

for i in range (n*n):
	f.write(str(random.randint(0, 100)))
	f.write(" ")
f.close()
