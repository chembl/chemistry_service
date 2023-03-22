```bash
mkdir etransafe_chemistry_services
cd etransafe_chemistry_services
conda create --name eTransafeChemistryServicesEnv
conda activate eTransafeChemistryServicesEnv
conda install -c anaconda django
conda install -c conda-forge djangorestframework
conda install -c conda-forge drf-yasg
conda install -c conda-forge chembl_structure_pipeline

django-admin startproject chemistry_services_project .

python manage.py startapp chemistry_services_app
python manage.py migrate
python manage.py createsuperuser 
```

