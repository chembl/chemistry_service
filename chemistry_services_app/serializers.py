from rest_framework import serializers
from django.core.serializers.json import Serializer


class QueryResultSerializer(serializers.Serializer):
    result = serializers.CharField()


class N2SQueryResultSerializer(Serializer):

    def get_dump_object(self, obj):

        if obj is None:
            return {}

        result_object = {
            'canonical_smiles': obj["canonical_smiles"],
            'standard_inchi': obj["standard_inchi"],
            'standard_inchi_key': obj["standard_inchi_key"]
        }

        return result_object


class SmilesRequestSerializer(serializers.Serializer):
    smiles = serializers.CharField()


class PAStandardizeRequestSerializer(serializers.Serializer):
    pa_type = serializers.CharField(required=False)
    compound = serializers.CharField()

