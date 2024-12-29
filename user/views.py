from django.shortcuts import render
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
            result = calculate_results(form.cleaned_data) 
            return render(request, 'calculation_results.html', {'result': result})
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