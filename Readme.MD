# eTransafe Chemistry Services

Django based web service wrapper for chembl_structure_pipeline(which is a set of ChEMBL protocols used to standardise and salt strip molecules. First used in ChEMBL 26.)

More details about **chembl_structure pipeline library**: 
- [source_code](https://github.com/chembl/ChEMBL_Structure_Pipeline.git)
- [wiki](https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki) 

## DATABASE SETTINGS
A chembl database required to run chemistry services. Update 'chemistry_services/settings/local.py' for database settings.

## Local Installation

With docker:
```bash
docker-compose up
```

## Usage
1) First run application with docker-compose up. For development environment docker-compose -f docker-compose-dev.yml
   
2) Visit the api schema from: [http://localhost:8000/v1/api](http://localhost:8000/api/)
