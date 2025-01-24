import json
from channels.generic.websocket import WebsocketConsumer
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from io import BytesIO
import base64

class FiltrationConsumer(WebsocketConsumer):
    def connect(self):
        self.accept()

    def disconnect(self, close_code):
        pass

    def receive(self, text_data):
        data = json.loads(text_data)
        self.filtration(data)

    def filtration(self, data):
        # Extract parameters from data if needed, or use defaults
        k0 = 1e-13
        myu0 = 1e-2
        P0 = 5e+5
        bet_y = 3e-10

        tMax = 3600

        lv = 0
        lp = 500
        betta = 1
        alf = 1

        bet = 1 + betta

        tau = 10
        h = 0.05
        L = 40
        T = int(tMax / tau) + 1
        N = int(L / h) + 1

        P = np.zeros((T, N))
        P[:, 0] = P0
        P[0, :] = 0
        P[:, -1] = 0
        U = np.zeros((T, N))
        kap = k0 / (myu0 * bet_y)

        for i in range(N):
            P[1, i] = P[0, i]

        kv = lv / ((tau ** bet) * gamma(3 - bet))
        kp = lp / ((tau ** alf) * (h ** 2) * gamma(2 - alf))

        kv1 = lv / ((tau ** betta) * gamma(2 - betta))  # velocity
        kp2 = lp / ((tau ** alf) * gamma(2 - alf))  # velocity

        A = kap / h ** 2 + kap * kp
        B = 1 / tau + kv + 2 * kap / h ** 2 + 2 * kap * kp
        C = kap / h ** 2 + kap * kp
        alpha = np.zeros(N)
        beta = np.zeros(N)
        beta[1] = P0

        for j in range(1, T - 1):
            for i in range(1, N - 1):
                sv = 0
                for k in range(1, j):
                    sv += (P[k + 1, i] - 2 * P[k, i] + P[k - 1, i]) * ((j - k + 1) ** (2 - bet) - (j - k) ** (2 - bet))

                sp1 = 0
                sp2 = 0
                sp3 = 0
                st = 0
                for k in range(j):
                    sp1 += (P[k + 1, i + 1] - P[k, i + 1]) * ((j - k + 1) ** (1 - alf) - (j - k) ** (1 - alf))
                    sp2 += (P[k + 1, i] - P[k, i]) * ((j - k + 1) ** (1 - alf) - (j - k) ** (1 - alf))
                    sp3 += (P[k + 1, i - 1] - P[k, i - 1]) * ((j - k + 1) ** (1 - alf) - (j - k) ** (1 - alf))
                    st += (U[k + 1, i] - U[k, i]) * ((j - k + 1) ** (1 - betta) - (j - k) ** (1 - betta))  # velocity

            for i in range(1, N - 1):
                F = (P[j, i] / tau - kv * sv + 2 * kv * P[j, i] - kv * P[j - 1, i] +
                     kp * kap * sp1 - kap * kp * P[j, i + 1] - 2 * kap * kp * sp2 +
                     2 * kap * kp * P[j, i] + kap * kp * sp3 - kap * kp * P[j, i - 1])

                alpha[i + 1] = C / (B - A * alpha[i])
                beta[i + 1] = (F + A * beta[i]) / (B - A * alpha[i])

            P[j + 1, -1] = 0
            for i in range(N - 2, -1, -1):
                P[j + 1, i] = alpha[i + 1] * P[j + 1, i + 1] + beta[i + 1]

            for i in range(N - 1):
                U[j + 1, i] = abs((-k0 * (P[j + 1, i + 1] - P[j + 1, i] +
                                             kp2 * (sp1 + P[j + 1, i + 1] - P[j, i + 1] - sp2 -
                                             P[j + 1, i] + P[j, i])) / (myu0 * h * (1 + kv1)) -
                                             (kv1 * st - kv1 * U[j, i]) / (1 + kv1)))

            U[j + 1, -1] = U[j + 1, -2]
            
        self.send_chart(P, T, N, h)

    def send_chart(self, P, T, N, h):
        plt.plot(np.arange(N) * h, P[T - 2, :])
        buf = BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
        buf.close()

        self.send(text_data=json.dumps({
            'chart': image_base64
        }))
