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
    k0 = float(data.get('permeability', 1e-13))
    myu0 = float(data.get('fluid_viscosity', 1e-2))
    P0 = float(data.get('pressure', 5e+5))
    bet_y = float(data.get('elastic_capacity_coefficient', 3e-10))
    tMax = float(data.get('maximum_time', 3600))
    lv = float(data.get('relaxation_time1', 1000))
    lp = float(data.get('relaxation_time2', 500))
    betta = float(data.get('fractional_derivative_time1', 0.9))
    alf = float(data.get('fractional_derivative_time2', 0.7))
    bet = 1 + betta
    tau = float(data.get('grid_step_direction1', 10))
    h = float(data.get('grid_step_direction2', 0.05))
    L = float(data.get('distance', 40))
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
                st += (U[k + 1, i] - U[k, i]) * ((j - k + 1) ** (1 - betta) - (j - k) ** (1 - betta))

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

    x = np.arange(N) * h
    y = P[T - 2, :]

    return {'x': x.tolist(), 'y': y.tolist()}


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