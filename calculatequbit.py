%matplotlib inline
import qiskit
from qiskit import IBMQ
from qiskit import Aer
from qiskit import QuantumCircuit, ClassicalRegister,QuantumRegister, QiskitError
from qiskit.quantum_info.operators import Operator
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
# |t, x, g, c, a>, with (t+x) having n qubits (index+pattern), having (n-1) qubits and c having 2 qubits
# Also, ancilla qubits (a) as support for mct gate

genome_file = open("yeast_chr1.txt", "r")
seq_x = genome_file.read()
print(seq_x)
genome_file.close()
seq_x = seq_x[0:32]
seq_y = "GAT"
Q_t = ceil(log2(len(seq_x)))
Q_x = len(seq_y)
Q_g = Q_t + Q_x - 1
Q_c = 2
Q_anc = 1
total_qubits = Q_t + Q_x + Q_g + Q_c + Q_anc
print(total_qubits)
## Initialization of IBM QX
#IBMQ.enable_account(’INSERT TOKEN HERE’)
#provider = IBMQ.get_provider()
