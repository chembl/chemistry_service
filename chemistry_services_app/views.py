from chembl_structure_pipeline import standardizer, checker
from drf_spectacular.utils import extend_schema

from drf_spectacular.openapi import OpenApiParameter, OpenApiTypes, OpenApiResponse
from rest_framework import status
from rest_framework.decorators import renderer_classes, action
from rest_framework.renderers import StaticHTMLRenderer
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework.generics import GenericAPIView
from django.http.response import HttpResponse
import json
from rdkit import Chem

from .chemistry_utils import convert_to_mol, get_standard_smiles, get_inchi_key, get_truncated_inchi_key
from .serializers import QueryResultSerializer, SmilesRequestSerializer, \
    PAStandardizeRequestSerializer, N2SQueryResultSerializer
from chemistry_services_app.db_utils import *

name_parameter = OpenApiParameter('name', OpenApiTypes.STR, OpenApiParameter.QUERY,
                                  description="Name for structure search")
inchi_key_parameter = OpenApiParameter('inchi_key', OpenApiTypes.STR, OpenApiParameter.QUERY,
                                       description="Inchi key for Structure to Name search")
smiles_parameter = OpenApiParameter('smiles', OpenApiTypes.STR, OpenApiParameter.QUERY, style="form",
                                    description="SMILES string of the compound")
pa_type_parameter = OpenApiParameter('pa_type', OpenApiTypes.STR, OpenApiParameter.QUERY, style="form",
                                     description="Accepted values: "
                                                 "preclinical, clinical. Default value is preclinical "
                                                 "if pa_type is empty")

compound_parameter = OpenApiParameter('compound', OpenApiTypes.STR, OpenApiParameter.QUERY, style="form",
                                      description="Compound should be either a name(such as paracetamol) when "
                                                  "pa_type is clinical, or a SMILES string when pa_type is preclinical")

invalid_request_message = {'Invalid Request': 'Request needs to have smiles parameter, check documentation.'}
invalid_compound_parameter_message = {'Invalid Request': 'Request needs to have compound parameter, '
                                                         'check documentation.'}
invalid_pa_type_parameter_message = {'Invalid Request': 'Request needs to have pa_type parameter, check documentation.'}
invalid_smiles_message = {'Invalid smiles': 'Request needs to have valid smiles.'}
invalid_n2s_request_message = {'Invalid Request': 'Request needs to have name parameter, check documentation.'}
invalid_s2n_request_message = {'Invalid Request': 'Request needs to have inchi_key parameter, check documentation.'}

request_is_valid_but_no_result_found_message = {'Empty response': 'Request is valid. But no result found for query.'}


class N2SView(APIView):


    @extend_schema(
        description='Tries to perform a partial search for a given name. '
                    'Name should be provided in `name` query parameter. '
                    'Returns a list of dictionaries for structures found for provided `name`. '
                    'Each item in the list contains (name, inchikey and smiles)',
        operation_id='Name to Structure',
        parameters=[name_parameter],

        responses={
            200: OpenApiResponse(description='Success'),
            404: OpenApiResponse(description='Invalid request/Resource not found'),
        }
    )
    def get(self, request, **kwargs):

        name = request.GET['name']
        if name is None:
            return Response(invalid_n2s_request_message, status=status.HTTP_404_NOT_FOUND)

        name = name.replace('-', '*')

        if '*' not in name:
            name += '*'

        words = name.split()

        search_term = ''
        for word in words:
            term_part = "+" + word
            search_term += term_part

        name = search_term

        rows = check_molecule_dictionary_full_text_search(name)

        if rows is not None and len(rows) == 1 and rows[0][0] is None:
            error_message = rows[0][5]
            return Response({'Empty response': error_message}, status=status.HTTP_200_OK)

        if rows is None or len(rows) == 0:
            rows = check_molecule_synonyms_full_text_search(name)

        if rows is not None and len(rows) == 1 and rows[0][0] is None:
            error_message = rows[0][5]
            return Response({'Empty response': error_message}, status=status.HTTP_200_OK)

        # If rows still empty return no content response
        if rows is None or len(rows) == 0:
            return Response(request_is_valid_but_no_result_found_message, status=status.HTTP_200_OK)

        # result_data = parse_n2s_result_rows_for_pa_standardize_query(rows)
        result_data = parse_n2s_result_rows_for_n2s_query(rows)
        result_json = json.dumps(result_data, default=lambda x: x.__dict__)

        return HttpResponse(result_json, status=status.HTTP_200_OK, content_type="application/json")
        # return Response(N2SQueryResultSerializer(result_data).data)


class S2NView(APIView):


    @extend_schema(
        description='Tries to perform a search for a structure(as inchi key). '
                    'Inchi key should be provided in `inchi_key` query parameter. '
                    'Returns a list of dictionaries for structures found for provided `inchi_key`. '
                    'Each item in the list contains (canonical_smiles, standard_inchi, standard_inchi_key and max_phase)',
        operation_id='Structure to Name',
        parameters=[inchi_key_parameter],

        responses={
            200: OpenApiResponse(description='Success'),
            404: OpenApiResponse(description='Invalid request/Resource not found'),
        }
    )
    def get(self, request, **kwargs):

        inchi_key = request.GET['inchi_key']
        if inchi_key is None:
            return Response(invalid_s2n_request_message, status=status.HTTP_404_NOT_FOUND)

        rows = structure_to_pref_name(inchi_key)

        if rows is not None and len(rows) == 1 and rows[0][0] is None:
            error_message = rows[0][5]
            return Response({'Empty response': error_message}, status=status.HTTP_200_OK)

        if rows is None or len(rows) == 0:
            rows = structure_to_synonym(inchi_key)

        if rows is not None and len(rows) == 1 and rows[0][0] is None:
            error_message = rows[0][5]
            return Response({'Empty response': error_message}, status=status.HTTP_200_OK)

        # If rows still empty return no content response
        if rows is None or len(rows) == 0:
            return Response(request_is_valid_but_no_result_found_message, status=status.HTTP_200_OK)

        result_data = parse_s2n_result_rows_for_s2n_query(rows)
        result_json = json.dumps(result_data, default=lambda x: x.__dict__)

        return HttpResponse(result_json, status=status.HTTP_200_OK, content_type="application/json")


class AliveView(APIView):
    @renderer_classes([StaticHTMLRenderer])
    def get(self, request, **kwargs):
        data = 'Chemistry Services alive'
        response = HttpResponse(data, content_type='text/plain; charset=UTF-8')
        return response


class ReadyView(GenericAPIView):
    @renderer_classes([StaticHTMLRenderer])
    def get(self, request, **kwargs):
        data = 'Chemistry Services ready'
        response = HttpResponse(data, content_type='text/plain; charset=UTF-8')
        return response


class StandardizeCompoundView(APIView):


    @extend_schema(

        description='Standardizes a given structure according to a set of predefined rules. '
                    'Compound should be provided in `smiles` query parameter. '
                    'Returns a tuple of InChIKey and standardized smiles. ',
        operation_id='Standardize compound',
        parameters=[smiles_parameter],
        request=SmilesRequestSerializer,
        responses={
            200: OpenApiResponse(response=QueryResultSerializer, description='Success'),
            404: OpenApiResponse(description='Invalid request/Resource not found'),
        },
    )
    @action(detail=False, methods=['post'])
    def post(self, request, **kwargs):
        serializer = SmilesRequestSerializer(request.POST)
        smiles = serializer.data.get('smiles', None)
        if smiles is None:
            return Response(invalid_request_message, status=status.HTTP_404_NOT_FOUND)
        mol_block = convert_to_mol(smiles)
        if mol_block is None:
            return Response(invalid_smiles_message, status=status.HTTP_400_BAD_REQUEST)

        std_mol_block = standardizer.standardize_molblock(mol_block)

        standardized_smiles = get_standard_smiles(std_mol_block)
        inchi_key = get_inchi_key(smiles)

        result = inchi_key, standardized_smiles

        result_data = {"result": result}
        return Response(QueryResultSerializer(result_data).data)


class PAStandardizeView(APIView):


    @extend_schema(
        description='Single endpoint to standardise molecules and to retrieve structures from names'
                    ' for the Primitive Adapters. ',
        operation_id='Standardize Compound For Primitive Adaptors',
        parameters=[pa_type_parameter, compound_parameter],
        request=PAStandardizeRequestSerializer,
        responses={
            200: OpenApiResponse(response=QueryResultSerializer, description='Success'),
            404: OpenApiResponse(description='Invalid request/Resource not found'),
        },
    )
    def post(self, request, **kwargs):

        serializer = PAStandardizeRequestSerializer(request.POST)
        pa_type = serializer.data.get('pa_type', None)

        if pa_type is None:
            pa_type = 'preclinical'

        compound = request.data.get('compound', None)

        if compound is None:
            return Response(invalid_compound_parameter_message, status=status.HTTP_404_NOT_FOUND)

        if pa_type == 'preclinical':
            try:
                mol = Chem.MolFromSmiles(compound)
                mol = standardizer.standardize_mol(mol)
                mol, _ = standardizer.get_parent_mol(mol)
                standardized_smiles = Chem.MolToSmiles(mol)
                inchi = Chem.MolToInchi(mol)
                inchikey = Chem.InchiToInchiKey(inchi)
                truncated_inchi_key = get_truncated_inchi_key(inchikey)
                result_list = []
                result_list.append(PAStandardizePreClinicalResult(standardized_smiles, inchikey, truncated_inchi_key))
                result_data = {"result": result_list}
                result_json = json.dumps(result_data, default=lambda x: x.__dict__)
            except:
                return Response({'Empty response': 'Error parsing the SMILES, please check your structure'},
                                status=status.HTTP_200_OK)
            return HttpResponse(result_json, status=status.HTTP_200_OK, content_type="application/json")
        elif pa_type == 'clinical':

            rows = check_molecule_dictionary(compound)

            if rows is not None and len(rows) == 1 and rows[0][0] is None:
                error_message = rows[0][5]
                return Response({'Empty response': error_message}, status=status.HTTP_200_OK)

            if rows is None or len(rows) == 0:
                rows = check_molecule_synonyms(compound)

            if rows is not None and len(rows) == 1 and rows[0][0] is None:
                error_message = rows[0][5]
                return Response({'Empty response': error_message}, status=status.HTTP_200_OK)

            # If rows still empty return no content response
            if rows is None or len(rows) == 0:
                return Response(request_is_valid_but_no_result_found_message, status=status.HTTP_200_OK)

            result_data = parse_n2s_result_rows_for_n2s_query(rows)
            result_json = json.dumps(result_data, default=lambda x: x.__dict__)

            return HttpResponse(result_json, status=status.HTTP_200_OK, content_type="application/json")
        else:
            return Response(invalid_pa_type_parameter_message, status.HTTP_404_NOT_FOUND)


class GetParentCompoundView(APIView):


    @extend_schema(
        description='Generates a parent structure based on a set of rules and '
                    'defined lists of salts and solvents.'
                    'Compound should be provided in `smiles` query parameter.  '
                    'Returns the standard smiles of parent structure. ',
        operation_id='Get parent of a given compound',
        parameters=[smiles_parameter],
        request=SmilesRequestSerializer,
        responses={
            200: OpenApiResponse(response=QueryResultSerializer, description='Success'),
            404: OpenApiResponse(description='Invalid request/Resource not found'),
        }
    )
    def post(self, request, **kwargs):
        serializer = SmilesRequestSerializer(request.POST)
        smiles = serializer.data.get('smiles', None)

        if smiles is None:
            return Response(invalid_request_message, status=status.HTTP_404_NOT_FOUND)
        mol_block = convert_to_mol(smiles)
        if mol_block is None:
            return Response(invalid_smiles_message, status=status.HTTP_400_BAD_REQUEST)

        std_mol_block, _ = standardizer.get_parent_molblock(mol_block)
        std_smiles = get_standard_smiles(std_mol_block)

        result_data = {"result": std_smiles}
        return Response(QueryResultSerializer(result_data).data)


class CheckCompoundView(APIView):


    @extend_schema(
        description='Assesses the quality of a structure. '
                    'It highlights specific features or issues in the structure that may need to be revised.'
                    'Together with the description of the issue, '
                    'the checker process returns a penalty score (between 0-9) which reflects '
                    'the seriousness of the issue (the higher the score, the more critical is the issue).'
                    'Compound should be provided in `smiles` query parameter.  '
                    'Returns a list of issues found by checker process.',

        operation_id='Check compound',
        parameters=[smiles_parameter],
        request=SmilesRequestSerializer,
        responses={
            200: OpenApiResponse(response=QueryResultSerializer, description='Success'),
            400: OpenApiResponse(description='Invalid smiles'),
            404: OpenApiResponse(description='Invalid request/Resource not found'),
        }
    )
    def post(self, request, **kwargs):
        serializer = SmilesRequestSerializer(request.POST)
        smiles = serializer.data.get('smiles', None)

        if smiles is None:
            return Response(invalid_request_message, status=status.HTTP_404_NOT_FOUND)
        mol_block = convert_to_mol(smiles)
        if mol_block is None:
            return Response(invalid_smiles_message, status=status.HTTP_400_BAD_REQUEST)

        issues = checker.check_molblock(mol_block)
        result_data = [{"result": issues}]
        return Response(QueryResultSerializer(result_data, many=True).data)
