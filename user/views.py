import json
import math
from django.http import JsonResponse
from django.shortcuts import render
from scipy.special import gamma
import numpy as np
from django.utils.timezone import now, timedelta
from .models import PageAccess
from .forms import CalculationForm1, CalculationForm2, CalculationForm3

def code(request):
    return render(request,"code.html")
def home(request):
    last_7_days = [now().date() - timedelta(days=i) for i in range(7)]
    access_counts = [
        PageAccess.objects.filter(timestamp__date=day).count() for day in last_7_days
    ]
    
    data = {
        "labels": [day.strftime("%Y-%m-%d") for day in last_7_days],
        "data": access_counts
    }
    return render(request, "index.html",{"chart_data": json.dumps(data)})

def calculate_results(form_data):
    return form_data
def calculation_form(request, form_id):
    if form_id == 1:
        form = CalculationForm1
    elif form_id == 2:
        form = CalculationForm2
    elif form_id == 3:
        form = CalculationForm3
    else:
        form = None
    if request.method == "POST": 
        form = form(request.POST) 
        if form.is_valid(): 
            data = calculate_results(form.cleaned_data) 
            return render(request, 'calculation_results'+ str(form_id) +'.html', {'data': data, 'form_id': form_id})
    context = {'form': form}
    return render(request, 'calculation_form.html', context)
def program(request):
    return render(request,'program.html')
def about(request):
    return render(request,'about.html')
def thesis(request):
    return render(request,'thesis.html')
def scopus(request):
    return render(request,'scopus.html')

def function1(data):
    k0 = float(data.get('permeability', 1e-13))
    myu0 = float(data.get('fluid_viscosity', 1e-2))
    P0 = float(data.get('pressure', 5e+5))
    bet_y = float(data.get('elastic_capacity_coefficient', 3e-10))
    tMax = float(data.get('maximum_time', 3600))
    lv = float(data.get('relaxation_time1', 0))
    lp = float(data.get('relaxation_time2', 500))
    alf = float(data.get('fractional_derivative_time1', 0.9))
    betta = float(data.get('fractional_derivative_time2', 0.7))
    tau = float(data.get('grid_step_direction1', 10))
    h = float(data.get('grid_step_direction2', 0.05))
    L = float(data.get('distance', 40))
    
    bet = 1 + betta
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

    kv1 = lv / ((tau ** betta) * gamma(2 - betta))
    kp2 = lp / ((tau ** alf) * gamma(2 - alf))

    A = kap / h ** 2 + kap * kp
    B = 1 / tau + kv + 2 * kap / h ** 2 + 2 * kap * kp
    C = kap / h ** 2 + kap * kp
    alpha = np.zeros(N)
    beta = np.zeros(N)
    beta[1] = P0

    for j in range(1, T - 1):
        for i in range(1, N - 1):
            sv = sum(
                (P[k + 1, i] - 2 * P[k, i] + P[k - 1, i]) * 
                ((j - k + 1) ** (2 - bet) - (j - k) ** (2 - bet)) 
                for k in range(1, j)
            )

        for i in range(1, N - 1):
            F = P[j, i] / tau - kv * sv + 2 * kv * P[j, i] - kv * P[j - 1, i]
            alpha[i + 1] = C / (B - A * alpha[i])
            beta[i + 1] = (F + A * beta[i]) / (B - A * alpha[i])

        P[j + 1, -1] = 0
        for i in range(N - 2, -1, -1):
            P[j + 1, i] = alpha[i + 1] * P[j + 1, i + 1] + beta[i + 1]

    x_values = np.arange(N) * h
    y_values = P[T - 2, :]

    # **NumPy massivlarini JSON ga mos formatga aylantirish**
    
    return {
        "x": list(x_values),  # `.tolist()` yoki `list()` orqali aylantirish
        "y": list(y_values)
    }

def function2(data):
    # Formadan kelgan parametrlarni o'qish
    k0 = float(data.get('permeability', 1e-14))
    myu0 = float(data.get('fluid_viscosity', 1.768e-3))
    P0 = float(data.get('pressure', 1e5))
    bet_y = float(data.get('elastic_capacity_coefficient', 3e-8))
    lv = float(data.get('relaxation_time1', 100))
    lp = float(data.get('relaxation_time2', 500))
    tau = float(data.get('grid_step_direction1', 10))
    h = float(data.get('grid_step_direction2', 0.05))
    alf = float(data.get('fractional_derivative_time1', 0.9))
    betta = float(data.get('fractional_derivative_time2', 0.7))
    L = float(data.get('distance', 5))
    tMax = float(data.get('maximum_time', 1800))
    g0 = float(data.get('gradient_pressure_limit', 1e5))

    bet = 1 + betta
    T = int(tMax / tau + 1)
    N = int(L / h + 1)

    dx = np.zeros((T, N))
    P = np.zeros((T, N))
    U = np.zeros((T, N))
    vaqt = np.zeros(N)

    kap = k0 / (myu0 * bet_y)

    P[:, 0] = P0
    P[0, :] = 0
    P[:, N - 1] = 0
    U[:, :] = 0

    for i in range(N):
        P[1, i] = P[0, i]

    kv = lv / (tau**bet * gamma(3 - bet))
    kp = lp / (tau**alf * h**2 * gamma(2 - alf))

    kv1 = lv / (tau**betta * gamma(2 - betta))
    kp2 = lp / (tau**alf * gamma(2 - alf))

    A = kap / h**2 + kap * kp
    B = 1 / tau + kv + 2 * kap / h**2 + 2 * kap * kp
    C = kap / h**2 + kap * kp
    alpha = np.zeros(N)
    beta = np.zeros(N)
    alpha[1] = 0
    beta[1] = P0

    for j in range(1, T - 1):
        for i in range(1, N - 1):
            sv = 0
            for k in range(1, j):
                sv += (P[k + 1, i] - 2 * P[k, i] + P[k - 1, i]) * ((j - k + 1)**(2 - bet) - (j - k)**(2 - bet))

            sp1 = 0
            sp2 = 0
            sp3 = 0
            st = 0
            for k in range(j):
                sp1 += (P[k + 1, i + 1] - P[k, i + 1]) * ((j - k + 1)**(1 - alf) - (j - k)**(1 - alf))
                sp2 += (P[k + 1, i] - P[k, i]) * ((j - k + 1)**(1 - alf) - (j - k)**(1 - alf))
                sp3 += (P[k + 1, i - 1] - P[k, i - 1]) * ((j - k + 1)**(1 - alf) - (j - k)**(1 - alf))
                st += (U[k + 1, i] - U[k, i]) * ((j - k + 1)**(1 - betta) - (j - k)**(1 - betta))

        for i in range(1, N - 1):
            F = (P[j, i] / tau - kv * sv + 2 * kv * P[j, i] - kv * P[j - 1, i] +
                 kp * kap * sp1 - kap * kp * P[j, i + 1] - 2 * kap * kp * sp2 +
                 2 * kap * kp * P[j, i] + kap * kp * sp3 - kap * kp * P[j, i - 1])

            alpha[i + 1] = C / (B - A * alpha[i])
            beta[i + 1] = (F + A * beta[i]) / (B - A * alpha[i])

        for i in range(N - 2, 0, -1):
            P[j + 1, i] = alpha[i + 1] * P[j + 1, i + 1] + beta[i + 1]

        for i in range(N - 1):
            dx[j + 1, i] = abs((P[j + 1, i + 1] - P[j + 1, i]) / h)

        for i in range(2, N - 1):
            if abs((P[j + 1, i] - P[j + 1, i - 1]) / h) > g0:
                if (P[j + 1, i] - P[j + 1, i - 1]) / h > 0:
                    U[j + 1, i] = abs((-k0 * (P[j + 1, i + 1] - P[j + 1, i] - g0 * h +
                                         kp2 * (sp1 + P[j + 1, i + 1] - P[j, i + 1] -
                                                sp2 - P[j + 1, i] + P[j, i])) /
                                        (myu0 * h * (1 + kv1))) - (kv1 * st - kv1 * U[j, i]) / (1 + kv1))
                    if j == T - 2:
                        print(i)
                else:
                    U[j + 1, i] = abs((-k0 * (g0 * h + (P[j + 1, i + 1] - P[j + 1, i]) +
                                         kp2 * (sp1 + P[j + 1, i + 1] - P[j, i + 1] -
                                                sp2 - P[j + 1, i] + P[j, i])) /
                                        (myu0 * h * (1 + kv1))) - (kv1 * st - kv1 * U[j, i]) / (1 + kv1))
            else:
                U[j + 1, i] = abs(-(kv1 * st - kv1 * U[j, i]) / (1 + kv1))

        U[j + 1, N - 1] = U[j + 1, N - 2]
        for i in range(N - 1):
            if dx[j + 1, i] >= g0 and vaqt[i] == 0:
                vaqt[i] = j * tau
    x_values = (np.arange(N) * h).tolist()  # X o‘qi (masofa)
    y_values = abs(U[T - 2, :]).tolist()  # Y o‘qi (modullangan U)

    result = {
        "x": x_values,
        "y": y_values
    }
    return result

def function3(data):
    k0 = float(data.get('permeability', 1e-13))
    myu0 = float(data.get('fluid_viscosity', 1e-2))
    Pc = float(data.get('constant_pressure', 150e5))
    Pk = float(data.get('initial_pressure', 200e5))
    bet_y = float(data.get('elastic_capacity_coefficient', 3e-10))
    lv = float(data.get('relaxation_time1', 1000))
    lp = float(data.get('relaxation_time2', 500))
    betta = float(data.get('fractional_derivative_time2', 0.7))
    alf = float(data.get('fractional_derivative_time1', 1))
    tau = float(data.get('grid_step_direction1', 100))
    h = float(data.get('grid_step_direction2', 0.05))
    L = float(data.get('distance', 40))
    tMax = float(data.get('maximum_time', 3600))
    rc = float(data.get('gradient_pressure_limit', 0.1))

    # Hisoblash
    bet = 1 + betta
    T = int(tMax / tau + 1)
    N = int(L / h + 1)
    Nc = math.floor(rc / h) + 1

    P = np.zeros((T, N))
    U = np.zeros((T, N))
    P[:, :Nc] = Pc
    P[0, Nc:] = Pk
    P[:, N - 1] = Pk
    U[:, :] = 0

    kap = k0 / (myu0 * bet_y)

    for i in range(3, N):
        P[1, i] = P[0, i]

    kv = lv / (tau**bet * gamma(3 - bet))
    kp = lp / (tau**alf * h**2 * gamma(2 - alf))

    kv1 = lv / (tau**betta * gamma(2 - betta))
    kp2 = lp / (tau**alf * gamma(2 - alf))

    for j in range(1, T - 1):
        for i in range(3, N - 1):
            sv = 0
            for k in range(1, j):
                sv += (P[k + 1, i] - 2 * P[k, i] + P[k - 1, i]) * ((j - k + 1)**(2 - bet) - (j - k)**(2 - bet))

            sp1 = 0
            sp2 = 0
            sp3 = 0
            st = 0
            for k in range(j):
                sp1 += (P[k + 1, i + 1] - P[k, i + 1]) * ((j - k + 1)**(1 - alf) - (j - k)**(1 - alf))
                sp2 += (P[k + 1, i] - P[k, i]) * ((j - k + 1)**(1 - alf) - (j - k)**(1 - alf))
                sp3 += (P[k + 1, i - 1] - P[k, i - 1]) * ((j - k + 1)**(1 - alf) - (j - k)**(1 - alf))
                st += (U[k + 1, i] - U[k, i]) * ((j - k + 1)**(1 - betta) - (j - k)**(1 - betta))

        alpha = np.zeros(N)
        beta = np.zeros(N)
        alpha[4] = 0
        beta[4] = Pc

        for i in range(4, N - 1):
            r1 = (2 * i * h - h) / 2
            r0 = (2 * i * h - 3 * h) / 2
            A = kap * (r0 / h**2 + kp * r0) / ((i - 1) * h)
            B = kap * ((r1 + r0) / h**2 + kp * (r1 + r0)) / ((i - 1) * h) + kv + 1 / tau
            C = kap * (r1 / h**2 + kp * r1) / ((i - 1) * h)
            F = P[j, i] / tau - kv * sv + 2 * kv * P[j, i] - kv * P[j - 1, i] + \
                (kap * kp * r1 * sp1) / ((i - 1) * h) - (kap * kp * r1 * P[j, i + 1]) / ((i - 1) * h) - \
                (kap * kp * (r0 + r1) * sp2) / ((i - 1) * h) + (kap * kp * (r1 + r0) * P[j, i]) / ((i - 1) * h) + \
                (kap * kp * r0 * sp3) / ((i - 1) * h) - (kap * kp * r0 * P[j, i - 1]) / ((i - 1) * h)

            alpha[i + 1] = C / (B - A * alpha[i])
            beta[i + 1] = (F + A * beta[i]) / (B - A * alpha[i])

        for i in range(N - 2, 3, -1):
            P[j + 1, i] = alpha[i + 1] * P[j + 1, i + 1] + beta[i + 1]

        for i in range(2, N - 1):
            U[j + 1, i] = abs((k0 * (P[j + 1, i + 1] - P[j + 1, i] + kp2 * (sp1 + P[j + 1, i + 1] - P[j, i + 1] - sp2 - P[j + 1, i] + P[j, i])) / (myu0 * h * (1 + kv1))) - (kv1 * st - kv1 * U[j, i]) / (1 + kv1))

        U[j + 1, N - 1] = U[j + 1, N - 2]

    x_values = (np.arange(N) * h).tolist()
    y_values = P[T - 2, :].tolist()

    result = {
        "x": x_values,
        "y": y_values
    }

    return result
def calculate(request, id):
    if request.method == 'POST':
        data = json.loads(request.POST.get('data'))
        if id == 1:
            result = function1(data)
        elif id == 2:
            result = function2(data)
        elif id == 3:
            result = function3(data)
        else:
            result = 'No valid function for this ID'
        return JsonResponse({'data': result})


def filtration_api(request):
    alf = float(request.GET.get('fractional_derivative_time1', 1))  
    betta = float(request.GET.get('fractional_derivative_time2', 1))  

    lambda_p = float(request.GET.get('relaxation_time2', 0))  
    lambda_v = float(request.GET.get('relaxation_time1', 0))

    result = filtration1(
        alf=alf,
        k0=float(request.GET.get('permeability', 1e-13)),
        myu0=float(request.GET.get('fluid_viscosity', 1e-2)),
        P0=float(request.GET.get('pressure', 5e+5)),
        bet_y=float(request.GET.get('elastic_capacity_coefficient', 3e-10)),
        lv=lambda_v,
        lp=lambda_p,
        tau=float(request.GET.get('grid_step_direction1', 10)),
        h=float(request.GET.get('grid_step_direction2', 0.05)),
        betta=betta
    )
    return JsonResponse({'result': result})

def filtration1(alf, tMax=3600, k0=1e-13, myu0=1e-2, P0=5e+5, bet_y=3e-10, lv=0, lp=500, betta=1, tau=10, h=0.05, L=40):
    bet = 1 + betta
    T = int(tMax / tau + 1)
    N = int(L / h + 1)
    P = np.zeros((T, N))
    U = np.zeros((T, N))
    P[:, 0] = P0
    kap = k0 / (myu0 * bet_y)
    kv = lv / ((tau ** bet) * gamma(3 - bet))
    kp = lp / ((tau ** alf) * (h ** 2) * gamma(2 - alf))
    kv1 = lv / ((tau ** betta) * gamma(2 - betta))
    kp2 = lp / ((tau ** alf) * gamma(2 - alf))
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
    
    return {"pressure": P[T - 2, :].tolist(), "velocity": (np.arange(N) * h).tolist()}