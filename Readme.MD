# eTransafe Chemistry Services
The chemistry service is a Django based web service wrapper for delivering compound identification functionalities. It is based on the ChEMBL database and the chembl_structure_pipeline (which is a set of ChEMBL protocols used to standardise and salt strip molecules - first used in ChEMBL 26). It was developed at the European Bioinformatic Institute as part of the eTransafe project. 

More details about **chembl_structure pipeline library**: 
- [source_code](https://github.com/chembl/ChEMBL_Structure_Pipeline.git)
- [wiki](https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki)
- [publication](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00456-1)

## Requirements
- Docker
- A MySQL instance of the ChEMBL database. ChEMBL can be freely [downloaded](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/).


## Local Installation Steps:
___
### Step 1:
Create .env file in the same folder with manage.py and put the line below:
```python
# (Note that you can set any value)
SECRET_KEY = 'n$i7dgz-qci*!fpecgbt9z$2onadzsl#za%!SDF_SADF_A3d'
```
### Step 2:
Update the contents of `chemistry_services_project/settings/local.py` with your ChEMBL mysql database credentials.
### Step 3:
To improve the performance of the chemistry service, it is advised to create full indexes on `pref_name` and `synonyms` in your ChEMBL database
```sql
"ALTER TABLE molecule_dictionary ADD FULLTEXT (pref_name); ALTER TABLE molecule_synonyms ADD FULLTEXT (synonyms);"
```
### Step 4:
```bash
docker-compose up --build
```
### Step 5:
- Visit the api schema from: [http://localhost:8000/v1/api](http://localhost:8000/v1/api)
___

## Usage

## Funding
This project has received funding from the Innovative Medicines Initiative 2 Joint Undertaking under grant
agreement No 777365 (“eTRANSAFE”)
