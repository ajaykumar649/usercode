####################################################
####################################################
## Author: Ajay Kumar                             ##
## Assistant Professor                            ##
## Sri Aurobindo College,University of Delhi      ##
####################################################
####################################################

import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

#fold = open('Poonam/175iri.txt', 'r') # Read from the input file
fnew1 = open('output1.txt','w') # Define output file to write
fnew2 = open('output2.txt','w') # Define output file to write

floats1 = []
floats2 = []


def func1(line):
    #fold = open('filename','r')
    #for line in filename:
        #floats1.extend([float(number) for number in line.split()])
        #Search for the line having DECOUPLING FACTOR
    if 'DECOUPLING FACTOR' in line:
        decp = line
        ddecp = decp.split(' ')
        #print ddecp
        dddecp = str(ddecp[-1:])
        ddddecp = dddecp.strip('\n')
        x = float(ddddecp[2:-4])
        fnew1.write(decp)
        return x
    ##y = x*2
            #print("x=", x)

        #print decp
    #Search for line having LEVEL NO and ENERGY

def func2(line):
    #fold = open('filename','r')
    #for line in filename:
    if 'LEVEL NO' in line:
        eng = line
        eeng = eng.split(' ')
        eeeng = str(eeng[-2:-1])
        y = float(eeeng[2:-2])
        fnew2.write(eng)
        return y
        #print("y", y)
        #print eng

#def caculate_E(jap):



#fold = open('/home/poonam/Desktop/forPoonam/175iri.txt', 'r') # Read from the input file
fold = open('1iri.txt', 'r') # Read from the input file
for line in fold:
    if (func1(line)):
        floats1.append(func1(line))
        #print x
    if (func2(line)):
        floats2.append(func2(line))
        #print func2(line)
#print ('Decoupling factor=',floats1)
print ('Energy=',floats2)

E0=0.0
h_sqr = 6.58*pow(10,-22)
g = 10.0
epsilon_o_k = E0-((h_sqr*h_sqr)/(2*g))
delta_k =1.0
spin = [2.5,3.5,4.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5]

floats3 =[]

for x in spin:
    #E_I_K = epsilon_o_k +  ((h_sqr*h_sqr) / (2 * g)) * (x *( x+1) ) + delta_k * floats1[0]*math.pow((-1),(x+1/2))*(x+1/2)
    E_I_K = epsilon_o_k +  ((h_sqr*h_sqr) / (2 * g)) * (x *( x+1) ) + delta_k * floats1[0]*cmath.exp(((x+1/2)**2)*cmath.log(-1))
    #cmath.exp(1.3 * cmath.log(-1.07))
    floats3.append(E_I_K)
    a = np.array(floats3)

#print ('coupling and energy', spin, floats3)

print a.real
print a.size
#print floats2
#print '%05f %05fi' % (n.real, n.imag)
print ('decoupling factor size is ',len(floats1))
print ('Energy read from diet file size is ',len(floats2))
print ('Energy read from diet file size is ',len(floats3))


plt.figure(1)
#plt.subplot(22)

plt.plot(spin, floats3)
#plt.yscale('linear')
plt.xlabel('I')
plt.ylabel('Energy')
plt.title('I vs Energy')
plt.grid(True)
plt.savefig("fig.png")

#fold.close()
#fnew1.close()
#fnew2.close()
