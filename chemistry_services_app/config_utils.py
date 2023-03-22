from django.conf import settings
import os

config_in_settings = settings.APPLICATION_CONFIG


def get_application_config_for_key(var_key):
    config_val = os.environ.get(var_key)
    if config_val is not None:
        print(var_key, config_val)
        return config_val
    elif var_key in config_in_settings:
        config_val = config_in_settings[var_key]
        print(var_key, config_val)
        return config_val
    else:
        return None






