import json
from django.http import JsonResponse
from django.shortcuts import render
from scipy.special import gamma
import numpy as np
from .forms import CalculationForm1, CalculationForm2, CalculationForm3

def home(request):
    return render(request, "index.html")

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
            return render(request, 'calculation_results.html', {'data': data, 'form_id': form_id})
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
    # O'zgaruvchilarni olish va boshlang'ich shartlar
    k0 = float(data.get('permeability', 1e-13))
    myu0 = float(data.get('fluid_viscosity', 1e-2))
    P0 = float(data.get('pressure', 5e+5))
    bet_y = float(data.get('elastic_capacity_coefficient', 3e-10))
    tMax = float(data.get('maximum_time', 3600))
    lv = float(data.get('relaxation_time1', 1000))
    lp = float(data.get('relaxation_time2', 500))
    betta = float(data.get('fractional_derivative_time1', 0.9))
    bet = 1 + betta
    tau = float(data.get('grid_step_direction1', 10))
    h = float(data.get('grid_step_direction2', 0.05))
    L = float(data.get('distance', 40))
    T = int(tMax / tau) + 1
    N = int(L / h) + 1

    alf_values = [1, 0.7, 0.5]
    results = {}

    kap = k0 / (myu0 * bet_y)

    # Hisoblashni takrorlash
    for alf in alf_values:
        P = np.zeros((T, N))
        P[:, 0] = P0
        U = np.zeros((T, N))

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

        # Asosiy hisoblash tsikli
        for j in range(1, T - 1):
            sv = np.zeros(N)
            sp1 = np.zeros(N)
            sp2 = np.zeros(N)
            sp3 = np.zeros(N)
            st = np.zeros(N)

            for k in range(1, j):
                delta = (j - k + 1) ** (2 - bet) - (j - k) ** (2 - bet)
                sv[1:-1] += (P[k + 1, 1:-1] - 2 * P[k, 1:-1] + P[k - 1, 1:-1]) * delta

            for k in range(j):
                delta_alf = (j - k + 1) ** (1 - alf) - (j - k) ** (1 - alf)
                delta_betta = (j - k + 1) ** (1 - betta) - (j - k) ** (1 - betta)

                sp1[1:-1] += (P[k + 1, 2:] - P[k, 2:]) * delta_alf
                sp2[1:-1] += (P[k + 1, 1:-1] - P[k, 1:-1]) * delta_alf
                sp3[1:-1] += (P[k + 1, :-2] - P[k, :-2]) * delta_alf
                st[1:-1] += (U[k + 1, 1:-1] - U[k, 1:-1]) * delta_betta

            # Tridiagonal sistemani hal qilish
            F = (P[j, 1:-1] / tau - kv * sv[1:-1] + 2 * kv * P[j, 1:-1] - kv * P[j - 1, 1:-1] +
                 kp * kap * sp1[1:-1] - kap * kp * P[j, 2:] - 2 * kap * kp * sp2[1:-1] +
                 2 * kap * kp * P[j, 1:-1] + kap * kp * sp3[1:-1] - kap * kp * P[j, :-2])

            alpha[2:] = C / (B - A * alpha[1:-1])
            beta[2:] = (F + A * beta[1:-1]) / (B - A * alpha[1:-1])

            P[j + 1, -1] = 0
            for i in range(N - 2, 0, -1):
                P[j + 1, i] = alpha[i + 1] * P[j + 1, i + 1] + beta[i + 1]

            # U ni hisoblash
            U[j + 1, :-1] = np.abs((-k0 * (P[j + 1, 1:] - P[j + 1, :-1] +
                                           kp2 * (sp1[1:] + P[j + 1, 1:] - P[j, 1:] - sp2[1:] -
                                                  P[j + 1, :-1] + P[j, :-1])) / (myu0 * h * (1 + kv1)) -
                                           (kv1 * st[1:] - kv1 * U[j, :-1]) / (1 + kv1)))
            U[j + 1, -1] = U[j + 1, -2]

        x = np.arange(N) * h
        y = P[T - 2, :]

        results[alf] = {'x': x.tolist(), 'y': y.tolist()}

    return results
def function2(data):
    return data

def function3(data):
    return data

def calculate(request, id):
    if request.method == 'POST':

        data = json.loads(request.POST.get('data'))
    
        print(type(data))
        if id == 1:
            result = function1(data)
        elif id == 2:
            result = function2(data)
        elif id == 3:
            result = function3(data)
        else:
            result = 'No valid function for this ID'
        return JsonResponse({'data': result})