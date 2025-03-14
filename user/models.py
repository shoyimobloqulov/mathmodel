from django.db import models

# Create your models here.
class PageAccess(models.Model):
    page_url = models.CharField(max_length=255)
    user_ip = models.GenericIPAddressField()
    timestamp = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.page_url} - {self.user_ip} - {self.timestamp}"
    
