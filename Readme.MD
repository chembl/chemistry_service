# eTransafe Chemistry Services
Django based web service wrapper for chembl_structure_pipeline(which is a set of ChEMBL protocols used to standardise and salt strip molecules. First used in ChEMBL 26.)

More details about **chembl_structure pipeline library**: 
- [source_code](https://github.com/chembl/ChEMBL_Structure_Pipeline.git)
- [wiki](https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki) 

## DATABASE SETTINGS
A chembl database required to run chemistry services. Update 'chemistry_services/settings/local.py' for database settings.

## Local Installation Steps:
___
### Step 1:
```python
# Create .env file in the same folder with manage.py and put the line below:
# (Note that you can set any value)
SECRET_KEY = 'n$i7dgz-qci*!fpecgbt9z$2onadzsl#za%!SDF_SADF_A3d'
```
### Step 2:
- Update the contents of `chemistry_services_project/settings/local.py` with your chembl mysql database credentials.
### Step 3:
```bash
docker-compose up --build
```
### Ste 4:
- Visit the api schema from: [http://localhost:8000/v1/api](http://localhost:8000/api/)
___