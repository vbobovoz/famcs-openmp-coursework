import matplotlib.pyplot as plt

K = [2, 4, 8, 16, 32, 64]
speedup_1 = [0.0190, 0.0204, 0.0509, 0.0737, 0.0440, 0.0721]
speedup_2 = [1.2874, 1.5341, 1.7547, 1.8521, 2.5462, 2.6315]
speedup_3 = [2.6163, 2.8241, 3.1003, 3.2294, 4.0831, 4.5938]

plt.figure(figsize=(10, 6))
plt.plot(K, speedup_1, marker='o', label='$N_1=9, N_2=11$')
plt.xticks(K)
plt.yticks(speedup_1)
plt.xlabel('Число итераций')
plt.ylabel('Ускорение')
plt.title('$N_1=9, N_2=11$')
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(K, speedup_2, marker='o', label='$N_1=99, N_2=101$')
plt.xticks(K)
plt.yticks(speedup_2)
plt.xlabel('Число итераций')
plt.ylabel('Ускорение')
plt.title('$N_1=99, N_2=101$')
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(K, speedup_3, marker='o', label='$N_1=999, N_2=1001$')
plt.xticks(K)
plt.yticks(speedup_3)
plt.xlabel('Число итераций')
plt.ylabel('Ускорение')
plt.title('$N_1=999, N_2=1001$')
plt.grid(True)
plt.show()