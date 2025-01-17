from django.urls import path
from . import views
urlpatterns = [
    path('', views.home, name='home'),
    path('topic/forms/<int:form_id>/', views.calculation_form, name='calculation_form'),
    path('program',views.program,name='program'),
    path('about',views.about,name="about"),
    path('thesis',views.thesis,name="thesis"),
    path('scopus',views.scopus,name='scopus')
]