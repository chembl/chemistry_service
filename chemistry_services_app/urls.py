from django.urls import include, path
from rest_framework import routers

from . import views

router = routers.DefaultRouter()

# Wire up our API using automatic URL routing.
# Additionally, we include login URLs for the browsable API.

urlpatterns = [
    path('', include(router.urls)),
    path('standardize_compound', views.StandardizeCompoundView.as_view(), name='standardize_compound'),
    path('get_parent_compound', views.GetParentCompoundView.as_view(), name='get_parent_compound'),
    path('check_compound', views.CheckCompoundView.as_view(), name='check_compound'),
    path('name_to_structure', views.N2SView.as_view(), name='name_to_structure'),
    path('structure_to_name', views.S2NView.as_view(), name='structure_to_name'),
    path('pa_standardize', views.PAStandardizeView.as_view(), name='pa_standardize'),
    path('alive', views.AliveView.as_view(), name='alive'),
    path('ready', views.ReadyView.as_view(), name='ready'),

]