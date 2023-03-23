# eTransafe Chemistry Services
Django based web service wrapper for chembl_structure_pipeline(which is a set of ChEMBL protocols used to standardise and salt strip molecules. First used in ChEMBL 26.)

More details about **chembl_structure pipeline library**: 
- [source_code](https://github.com/chembl/ChEMBL_Structure_Pipeline.git)
- [wiki](https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki) 

## REQUIREMENTS
- Docker
- A MySQL instance of the ChEMBL database. ChEMBL can be freely [downloaded](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/).


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
- To improve the performance of the chemistry service, it is advised to create full indexes on `pref_name` and `synonyms` in your ChEMBL database
```sql
"ALTER TABLE molecule_dictionary ADD FULLTEXT (pref_name); ALTER TABLE molecule_synonyms ADD FULLTEXT (synonyms);"
```
### Step 4:
```bash
docker-compose up --build
```
### Ste 5:
- Visit the api schema from: [http://localhost:8000/v1/api](http://localhost:8000/v1/api/)
___

## FUNDING
This project has received funding from the Innovative Medicines Initiative 2 Joint Undertaking under grant
agreement No 777365 (“eTRANSAFE”)
