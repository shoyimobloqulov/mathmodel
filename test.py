import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

# Dastlabki parametrlar
k0 = 1e-13
myu0 = 1e-2
P0 = 5e5
bet_y = 3e-10
tMax = 3600

lv = 0
lp = 500
betta = 1
alf = 1

bet = 1 + betta

tau = 100
h = 0.05
L = 40
T = int(tMax / tau) + 1
N = int(L / h) + 1

# Boshlang'ich va chegaraviy shartlar
P = np.zeros((T, N))
P[:, 0] = P0
U = np.zeros((T, N))
kap = k0 / (myu0 * bet_y)

for i in range(N):
    P[1, i] = P[0, i]

kv = lv / ((tau ** bet) * sp.gamma(3 - bet))
kp = lp / ((tau ** alf) * (h ** 2) * sp.gamma(2 - alf))
kv1 = lv / ((tau ** betta) * sp.gamma(2 - betta))
kp2 = lp / ((tau ** alf) * sp.gamma(2 - alf))

A = kap / h ** 2 + kap * kp
B = 1 / tau + kv + 2 * kap / h ** 2 + 2 * kap * kp
C = kap / h ** 2 + kap * kp

alpha = np.zeros(N)
beta = np.zeros(N)
beta[1] = P0

for j in range(1, T - 1):
    for i in range(1, N - 1):
        sv = sum((P[k + 1, i] - 2 * P[k, i] + P[k - 1, i]) * ((j - k + 1) ** (2 - bet) - (j - k) ** (2 - bet)) for k in range(1, j))
        sp1 = sum((P[k + 1, i + 1] - P[k, i + 1]) * ((j - k + 1) ** (1 - alf) - (j - k) ** (1 - alf)) for k in range(j))
        sp2 = sum((P[k + 1, i] - P[k, i]) * ((j - k + 1) ** (1 - alf) - (j - k) ** (1 - alf)) for k in range(j))
        sp3 = sum((P[k + 1, i - 1] - P[k, i - 1]) * ((j - k + 1) ** (1 - alf) - (j - k) ** (1 - alf)) for k in range(j))
        st = sum((U[k + 1, i] - U[k, i]) * ((j - k + 1) ** (1 - betta) - (j - k) ** (1 - betta)) for k in range(j))
    
    for i in range(1, N - 1):
        F = (P[j, i] / tau - kv * sv + 2 * kv * P[j, i] - kv * P[j - 1, i] + kp * kap * sp1 - kap * kp * P[j, i + 1] - 2 * kap * kp * sp2 + 2 * kap * kp * P[j, i] + kap * kp * sp3 - kap * kp * P[j, i - 1])
        alpha[i + 1] = C / (B - A * alpha[i])
        beta[i + 1] = (F + A * beta[i]) / (B - A * alpha[i])

    P[j + 1, N - 1] = 0
    for i in range(N - 2, 0, -1):
        P[j + 1, i] = alpha[i + 1] * P[j + 1, i + 1] + beta[i + 1]
    
    for i in range(N - 1):
        U[j + 1, i] = abs((-k0 * (P[j + 1, i + 1] - P[j + 1, i] + kp2 * (sp1 + P[j + 1, i + 1] - P[j, i + 1] - sp2 - P[j + 1, i] + P[j, i])) / (myu0 * h * (1 + kv1))) - (kv1 * st - kv1 * U[j, i]) / (1 + kv1))
    
    U[j + 1, N - 1] = U[j + 1, N - 2]

# Natijalarni chiqarish
plt.plot(np.arange(N) * h, P[T - 2, :], label='Bosim P')
plt.xlabel('Masofa')
plt.ylabel('Bosim')
plt.legend()
plt.grid()
plt.show()

print(np.arange(N) * h)
