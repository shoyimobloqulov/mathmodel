from django.urls import re_path
from .consumers import FiltrationConsumer

websocket_urlpatterns = [
    re_path(r'ws/filtration/$', FiltrationConsumer.as_asgi()),
]
