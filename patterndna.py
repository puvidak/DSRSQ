##Quantum Pattern Recognition
## Python code
### Original articles:
###
### (1) "Improving the Sequence Alignment Method by Quantum
Multi-Pattern Recognition"
### Konstantinos Prousalis & Nikos Konofaos
### Published in: SETN ’18 Proceedings of the 10th Hellenic
Conference on Artificial Intelligence, Article No. 50
### Patras, Greece, July 09 - 12, 2018
###
### (2) "Quantum Pattern Recognition with Probability of 100%"
### Rigui Zhou & Qiulin Ding
### Published in: International Journal of Theoretical Physics,
Springer
### Received: 3 August 2007, Accepted: 11 September 2007,
Published online: 4 October 2007
###
### (3) "Initializing the amplitude distribution of a quantum
state"
### Dan Ventura & Tony Martinez
### Revised 2nd November 1999
## Importing libraries
%matplotlib inline
import qiskit
from qiskit import IBMQ
from qiskit import Aer
from qiskit import QuantumCircuit, ClassicalRegister,
QuantumRegister, QiskitError
from qiskit.quantum_info.operators import Operator
101
C – Quantum Pattern Recognition
from qiskit.tools.visualization import circuit_drawer
from qiskit.tools.visualization import plot_histogram
from qiskit.tools.visualization import plot_state_city
from qiskit.providers.aer import noise
import random
from math import *
import math
import matplotlib
import numpy as np
## Initializing global variables
# Quantum register is organized like the following:
# |t, x, g, c, a>, with (t+x) having n qubits (index+pattern), g
having (n-1) qubits and c having 2 qubits
# Also, ancilla qubits (a) as support for mct gate
genome_file = open("HIVgenome.txt", "r")
seq_x = genome_file.read()
genome_file.close()
seq_x = seq_x[0:32]
seq_y = "GAT"
Q_t = ceil(log2(len(seq_x)))
Q_x = len(seq_y)
Q_g = Q_t + Q_x - 1
Q_c = 2
Q_anc = 1
total_qubits = Q_t + Q_x + Q_g + Q_c + Q_anc
## Initialization of IBM QX
IBMQ.enable_account(’INSERT TOKEN HERE’)
provider = IBMQ.get_provider()
# Pick an available backend
# If this isn’t available pick a backend whose name containes
’_qasm_simulator’ from the output above
backend = provider.get_backend(’ibmq_qasm_simulator’)
# Uncomment if you want to use local simulator
#backend= Aer.get_backend(’qasm_simulator’)
## Functions for recurrence dot matrix
def delta(x, y):
102
C – Quantum Pattern Recognition
return 0 if x == y else 1
def M(seq1, seq2, i, j, k):
return sum(delta(x, y) for x, y in zip(seq1[i : i+k],seq2[j :
j+k]))
def makeMatrix(seq1, seq2, k):
n = len(seq1)
m = len(seq2)
return [[M(seq1, seq2, i, j, k) for j in range(m - k + 1)]
for i in range(n - k + 1)]
def plotMatrix(M, t, seq1, seq2, nonblank = chr(0x25A0), blank =
’ ’):
print(’ |’ + seq2)
print(’-’ * (2 + len(seq2)))
for label, row in zip(seq1, M):
line = ’’.join(nonblank if s < t else blank for s in row)
print(label + ’|’ + line)
return
def convertMatrix(M):
for i in range(0, len(M)):
for j in range(0, len(M[i])):
if M[i][j] == 0:
M[i][j] = 1
elif M[i][j] == 1:
M[i][j] = 0
return M
def dotplot(seq1, seq2, k = 1, t = 1):
if len(seq1) > len(seq2):
raise Exception("Vertical sequence cannot be longer than
horizontal sequence!")
M = makeMatrix(seq1, seq2, k)
plotMatrix(M, t, seq1, seq2)
M = convertMatrix(M)
return M
def getAllDiagonalsFromMatrix(M):
D = np.array([])
d_size = -1
for i in range(0, len(M[0])):
d = np.diag(M, k=i)
if d_size == -1:
d_size = len(d)
D = d
103
C – Quantum Pattern Recognition
elif d_size > len(d):
z = np.zeros((1, (d_size-len(d))), dtype=int)
d = np.append(d, z)
D = np.vstack((D, d))
else:
D = np.vstack((D, d))
return D
def convertBinArrayToStr(array):
string = ""
for bin_digit in array:
if bin_digit == 0:
string = string + ’0’
elif bin_digit == 1:
string = string + ’1’
return string
## Functions for Quantum Pattern Recognition
def generateInitialState(qc, qr, dot_matrix):
D = getAllDiagonalsFromMatrix(dot_matrix)
m = len(D)
print("Size of Learning Set: {}".format(len(D)))
idx = 0
for d in D:
print("{}->{}: {}".format(idx, format(idx,
’0’+str(Q_t)+’b’), d))
idx = idx + 1
z_values = convertBinArrayToStr(np.zeros(Q_t+Q_x))
ancilla_qubits = []
for qi in range(0, Q_anc):
ancilla_qubits.append(qr[Q_t + Q_x + Q_g + Q_c + qi])
for p in range(m, 0, -1):
bin_diagonal = convertBinArrayToStr(D[len(D)-p])
index = format((len(D)-p), ’0’ + str(Q_t) + ’b’)
instance = index + bin_diagonal
#print("Instance #{}, z={}".format(p, instance))
for j in range(1, Q_t + Q_x + 1):
if z_values[j-1] != instance[j-1]:
#print("F_0 #{} Applied to circuit with ctrl={}
and target={}".format(j, Q_t+Q_x+Q_g+Q_c-1,
j-1))
qc.x(qr[Q_t+Q_x+Q_g+Q_c-1])
qc.cx(qr[Q_t+Q_x+Q_g+Q_c-1], qr[j-1])
qc.x(qr[Q_t+Q_x+Q_g+Q_c-1])
104
C – Quantum Pattern Recognition
z_values = instance
#print("F_0 Applied to circuit with ctrl={} and
target={}".format(Q_t+Q_x+Q_g+Q_c-1,
Q_t+Q_x+Q_g+Q_c-2))
qc.x(qr[Q_t+Q_x+Q_g+Q_c-1])
qc.cx(qr[Q_t+Q_x+Q_g+Q_c-1], qr[Q_t+Q_x+Q_g+Q_c-2])
qc.x(qr[Q_t+Q_x+Q_g+Q_c-1])
#print("S_{},{} Applied to circuit with ctrl={} and
target={}".format(1, p, Q_t+Q_x+Q_g+Q_c-2,
Q_t+Q_x+Q_g+Q_c-1))
theta = 2*np.arcsin(1/sqrt(p))
qc.cry(theta, qr[Q_t+Q_x+Q_g+Q_c-2],
qr[Q_t+Q_x+Q_g+Q_c-1])
if instance[0]==’0’ and instance[1]==’0’:
#print("A_00 #1 Applied to circuit with ctrl={},{} and
target={}".format(0, 1, Q_t+Q_x))
qc.x(qr[0])
qc.x(qr[1])
qc.mct([qr[0], qr[1]], qr[Q_t+Q_x], ancilla_qubits)
qc.x(qr[1])
qc.x(qr[0])
elif instance[0]==’0’ and instance[1]==’1’:
#print("A_01 #1 Applied to circuit with ctrl={},{} and
target={}".format(0, 1, Q_t+Q_x))
qc.x(qr[0])
qc.mct([qr[0], qr[1]], qr[Q_t+Q_x], ancilla_qubits)
qc.x(qr[0])
elif instance[0]==’1’ and instance[1]==’0’:
#print("A_10 #1 Applied to circuit with ctrl={},{} and
target={}".format(0, 1, Q_t+Q_x))
qc.x(qr[1])
qc.mct([qr[0], qr[1]], qr[Q_t+Q_x], ancilla_qubits)
qc.x(qr[1])
elif instance[0]==’1’ and instance[1]==’1’:
#print("A_11 #1 Applied to circuit with ctrl={},{} and
target={}".format(0, 1, Q_t+Q_x))
qc.mct([qr[0], qr[1]], qr[Q_t+Q_x], ancilla_qubits)
for k in range(3, Q_t+Q_x+1):
if instance[k-1]==’0’:
#print("A_01 #{} Applied to circuit with
ctrl={},{} and target={}".format(k-1, k-1,
Q_t+Q_x+k-3, Q_t+Q_x+k-2))
qc.x(qr[k-1])
105
C – Quantum Pattern Recognition
qc.mct([qr[k-1], qr[Q_t+Q_x+k-3]],
qr[Q_t+Q_x+k-2], ancilla_qubits)
qc.x(qr[k-1])
elif instance[k-1]==’1’:
#print("A_11 #{} Applied to circuit with
ctrl={},{} and target={}".format(k-1, k-1,
Q_t+Q_x+k-3, Q_t+Q_x+k-2))
qc.mct([qr[k-1], qr[Q_t+Q_x+k-3]],
qr[Q_t+Q_x+k-2], ancilla_qubits)
#print("F_1 Applied to circuit with ctrl={} and
target={}".format(Q_t+Q_x+Q_g-1, Q_t+Q_x+Q_g))
qc.cx(qr[Q_t+Q_x+Q_g-1], qr[Q_t+Q_x+Q_g])
for k in range(Q_t+Q_x, 2, -1):
if instance[k-1]==’0’:
#print("A_01 #{} Applied to circuit with
ctrl={},{} and target={}".format(k-1, k-1,
Q_t+Q_x+k-3, Q_t+Q_x+k-2))
qc.x(qr[k-1])
qc.mct([qr[k-1], qr[Q_t+Q_x+k-3]],
qr[Q_t+Q_x+k-2], ancilla_qubits)
qc.x(qr[k-1])
elif instance[k-1]==’1’:
#print("A_11 #{} Applied to circuit with
ctrl={},{} and target={}".format(k-1, k-1,
Q_t+Q_x+k-3, Q_t+Q_x+k-2))
qc.mct([qr[k-1], qr[Q_t+Q_x+k-3]],
qr[Q_t+Q_x+k-2], ancilla_qubits)
if instance[0]==’0’ and instance[1]==’0’:
#print("A_00 #1 Applied to circuit with ctrl={},{} and
target={}".format(0, 1, Q_t+Q_x))
qc.x(qr[0])
qc.x(qr[1])
qc.mct([qr[0], qr[1]], qr[Q_t+Q_x], ancilla_qubits)
qc.x(qr[1])
qc.x(qr[0])
elif instance[0]==’0’ and instance[1]==’1’:
#print("A_01 #1 Applied to circuit with ctrl={},{} and
target={}".format(0, 1, Q_t+Q_x))
qc.x(qr[0])
qc.mct([qr[0], qr[1]], qr[Q_t+Q_x], ancilla_qubits)
qc.x(qr[0])
elif instance[0]==’1’ and instance[1]==’0’:
#print("A_10 #1 Applied to circuit with ctrl={},{} and
target={}".format(0, 1, Q_t+Q_x))
qc.x(qr[1])
106
C – Quantum Pattern Recognition
qc.mct([qr[0], qr[1]], qr[Q_t+Q_x], ancilla_qubits)
qc.x(qr[1])
elif instance[0]==’1’ and instance[1]==’1’:
#print("A_11 #1 Applied to circuit with ctrl={},{} and
target={}".format(0, 1, Q_t+Q_x))
qc.mct([qr[0], qr[1]], qr[Q_t+Q_x], ancilla_qubits)
#print("F Applied to circuit at
qubit={}".format(Q_t+Q_x+Q_g+Q_c-1))
qc.x(qr[Q_t+Q_x+Q_g+Q_c-1])
return
def getIndices(mySet):
indices = []
for i in range(0, len(mySet)):
tmp = ""
for j in range(0, len(mySet[i])):
tmp = tmp + str(int(mySet[i][j]))
indices.append(int(tmp, 2))
return indices
def oracle(query_set):
I = np.identity(2**Q_x)
b_sum = np.zeros((2**Q_x, 2**Q_x))
indices = getIndices(query_set)
for i in indices:
vs = np.zeros((1, 2**Q_x))
for j in range(0, 2**Q_x):
if j == i:
vs[0][j] = 1
b_sum = b_sum + np.dot(np.conjugate(np.transpose(vs)), vs)
U = I - (1-1j)*b_sum
return U
def inversionAboutMean(dot_matrix):
I = np.identity(2**(Q_t+Q_x))
b_sum = np.zeros((2**(Q_t+Q_x), 1))
D = getAllDiagonalsFromMatrix(dot_matrix)
mySet = np.empty([len(D), Q_t+Q_x])
for i in range(0, len(D)):
bin_arr = np.concatenate((convertIntToBinArray(i, Q_t),
D[i]))
mySet[i] = bin_arr
indices = getIndices(mySet)
107
C – Quantum Pattern Recognition
for i in indices:
vs = np.zeros((2**(Q_t+Q_x), 1))
for j in range(0, 2**(Q_t+Q_x)):
if j == i:
vs[j][0] = 1
b_sum = b_sum + vs
phi_zero = (1/sqrt(len(D))) * b_sum
phi_mtrx = np.dot(phi_zero,
np.conjugate(np.transpose(phi_zero)))
U = (1 + 1j) * phi_mtrx - 1j * I
return U
def convertIntToBinArray(j, dim):
if not isinstance(j, int):
raise Exception("Number of bits must be an integer!")
elif (j == 0 or j == 1) and dim < 1:
raise Exception("More bits needed to convert j in
binary!")
elif j > 1 and dim <= log2(j):
raise Exception("More bits needed to convert j in
binary!")
bin_arr = np.array([], dtype=int)
j_bin = format(int(j), ’0’ + str(dim) + ’b’)
for k in j_bin:
if k == ’1’:
bin_arr = np.append(bin_arr, 1)
elif k == ’0’:
bin_arr = np.append(bin_arr, 0)
return bin_arr
def QPR(dot_matrix):
qr = qiskit.QuantumRegister(total_qubits)
cr = qiskit.ClassicalRegister(Q_t)
qc = qiskit.QuantumCircuit(qr, cr)
print("Total number of qubits: {}".format(total_qubits))
print("Size of t register: {}".format(Q_t))
print("Size of x register: {}".format(Q_x))
print("Size of g register: {}".format(Q_g))
print("Size of c register: {}".format(Q_c))
print("Number of ancilla qubits: {}".format(Q_anc))
# A query set is manually defined
query_set = np.array([[1,1,1],
108
C – Quantum Pattern Recognition
[0,1,1],
[1,1,0],
[1,0,1]])
O_mtrx = oracle(query_set)
U_phi_mtrx = inversionAboutMean(dot_matrix)
O = Operator(O_mtrx)
U_phi = Operator(U_phi_mtrx)
O_qubits = []
for qi in range(Q_x-1, -1, -1):
O_qubits.append(Q_t + qi)
U_phi_qubits = []
for qi in range(Q_t+Q_x-1, -1, -1):
U_phi_qubits.append(qi)
generateInitialState(qc, qr, dot_matrix)
#simulateStateVector(qc)
T = round((math.pi/4)*sqrt(len(dot_matrix[0])/len(query_set)))
it = 0
for i in range(0, T):
print("Grover Iteration #{}".format(it+1))
qc.unitary(O, O_qubits, label=’O’)
#simulateStateVector(qc)
qc.unitary(U_phi, U_phi_qubits, label=’U_phi’)
#simulateStateVector(qc)
it = it + 1
print("Grover’s algorithm had {} iterations.".format(int(it)))
finalGroverMeasurement(qc, qr, cr)
return qc
def simulateStateVector(qc):
result = qiskit.execute(qc,
backend=Aer.get_backend(’statevector_simulator’)).result()
state = result.get_statevector(qc)
print("Number of states in vector: {}".format(len(state)))
it = 0
for item in state:
bin_str = format(it, ’0’+str(total_qubits)+’b’)
bin_str_rev = bin_str[len(bin_str)::-1]
if (item.real != 0 or item.imag != 0):
print("{}->{}: {}".format(it,
bin_str_rev[Q_t:Q_t+Q_x], item))
109
C – Quantum Pattern Recognition
it = it + 1
return
# Final measurement
def finalGroverMeasurement(qc, qr, cr):
for qi in range(0, Q_t):
qc.measure(qr[qi], cr[qi])
return
## Main function
if __name__ == ’__main__’:
# Printing some data for testing
M = dotplot(seq_y, seq_x)
qc = QPR(M)
print("Circuit depth: {}".format(qc.depth()))
# Total number of gates
print("Number of gates: {}".format(len(qc.data)))
gate_num = 1
for item in qc.data:
qb_list = ""
for qb in item[1]:
qb_list = qb_list + str(qb.index) + ", "
qb_list = qb_list[:len(qb_list)-2]
print("#{}: {}, {}".format(gate_num, item[0].name,
qb_list))
gate_num = gate_num + 1
# Drawing circuit
#qc.draw()
# Showing histogram
# BE CAREFUL!
# Qiskit uses a LSB ordering, meaning the first qubit is all
the way to the right!
# For example, a state of |01> would mean the first qubit is
1 and the second qubit is 0!
sim = qiskit.execute(qc, backend=backend, shots=1024)
result = sim.result()
final=result.get_counts(qc)
print(final)
plot_histogram(final)
