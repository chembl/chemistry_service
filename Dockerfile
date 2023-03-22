FROM continuumio/miniconda3

WORKDIR /usr/src/app
COPY . /usr/src/app/

# Create the environment:
RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y curl
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "eTransafeEnv", "/bin/bash", "-c"]

RUN pip install drf-spectacular

# Make sure the environment is activated:
RUN echo "Make sure strandardizer, django and rdkit is installed:"
RUN python -c "import chembl_structure_pipeline"
RUN python -c "import django"
RUN python -c "import rest_framework"
RUN python -c "import rdkit"
RUN python -c "import drf_spectacular"
ENV PYTHONUNBUFFERED=1

EXPOSE 8000
ENTRYPOINT ["python", "manage.py", "runserver", "0.0.0.0:8000", "--settings=chemistry_services_project.settings.local"]
