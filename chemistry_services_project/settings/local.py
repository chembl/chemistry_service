from .base import *

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'db.sqlite3',
    },
    'chembl_mysql': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': os.environ.get('MYSQL_DB_NAME', 'chembl_28'),
        'USER': os.environ.get('MYSQL_USER', 'chembl_28'),
        'PASSWORD': os.environ.get('MYSQL_PASSWORD', 'chembl_28'),
        'HOST': os.environ.get('MYSQL_HOST', '#ENTER_MYSQL_HOSTNAME#'),
        'PORT': os.environ.get('MYSQL_PORT', '#ENTER_MYSQL_PORT#')
    }
}