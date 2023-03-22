from django.db import connections
from .chemistry_utils import get_truncated_inchi_key


def check_molecule_dictionary(query):
    with connections['chembl_mysql'].cursor() as cursor:

        rows_count = cursor.execute(
            "select distinct cs.canonical_smiles, cs.standard_inchi, cs.standard_inchi_key, md.max_phase, "
            "case "
            "  when cr.src_id in (8, 9, 12, 13, 36, 41, 42, 53) then 1 "
            "  else 2 "
            "  end as src_priority, "
            "case "
            " when md.structure_type = 'NONE' and cp.full_molformula is not null then 'Metal-containing compound' "
            " else md.molecule_type end molecule_type, "
            "md.pref_name "
            "from compound_records cr "
            "    join molecule_hierarchy mh on cr.molregno = mh.molregno "
            "    left join compound_structures cs on mh.parent_molregno = cs.molregno "
            "    left join molecule_dictionary md on mh.parent_molregno = md.molregno "
            "    left join compound_properties cp ON mh.parent_molregno = cp.molregno "
            "where md.pref_name is not null "
            "	and cr.src_id in (1, 8, 9, 12, 13, 36, 41, 42, 53) "
            "and md.pref_name = %s"
            "order by md.max_phase desc, src_priority limit 1", [query])

        if rows_count > 0:
            rows = cursor.fetchall()
        else:
            return None

    return rows


def check_molecule_synonyms(query):
    with connections['chembl_mysql'].cursor() as cursor:

        rows_count = cursor.execute("select cs.canonical_smiles, cs.standard_inchi, cs.standard_inchi_key, md.max_phase, "
                                    "case "
                                    "  when cr.src_id in (8, 9, 12, 13, 36, 41, 42, 53) then 1 "
                                    "  else 2 "
                                    "end as src_priority, "
                                    "case "
                                    " when md.structure_type = 'NONE' and cp.full_molformula is not null then 'Metal-containing compound' "
                                    " else md.molecule_type end molecule_type, "
                                    "upper(ms.synonyms) "
                                    "from molecule_synonyms ms "
                                    " join compound_records cr on ms.molregno = cr.molregno "
                                    " join molecule_hierarchy mh on cr.molregno = mh.molregno "
                                    " join compound_structures cs on mh.parent_molregno = cs.molregno "
                                    " join molecule_dictionary md on mh.parent_molregno = md.molregno "
                                    " left join compound_properties cp ON mh.parent_molregno = cp.molregno "
                                    "where cr.src_id in (1, 8, 9, 12, 13, 36, 41, 42, 53) "
                                    " and ms.synonyms = %s "
                                    "group by cs.canonical_smiles, cs.standard_inchi, cs.standard_inchi_key, md.max_phase, src_priority "
                                    "order by md.max_phase desc, src_priority limit 1", [query])

        if rows_count > 0:
            rows = cursor.fetchall()
        else:
            return None

    return rows


def check_molecule_dictionary_full_text_search(query):
    with connections['chembl_mysql'].cursor() as cursor:

        rows_count = cursor.execute(
            "select distinct cs.canonical_smiles, cs.standard_inchi, cs.standard_inchi_key, md.max_phase, "
            "case "
            "  when cr.src_id in (8, 9, 12, 13, 36, 41, 42, 53) then 1 "
            "  else 2 "
            "  end as src_priority, "
            "case "
            " when md.structure_type = 'NONE' and cp.full_molformula is not null then 'Metal-containing compound' "
            " else md.molecule_type end molecule_type, "
            "md.pref_name "
            "from compound_records cr "
            "    join molecule_hierarchy mh on cr.molregno = mh.molregno "
            "    left join compound_structures cs on mh.parent_molregno = cs.molregno "
            "    left join molecule_dictionary md on mh.parent_molregno = md.molregno "
            "    left join compound_properties cp ON mh.parent_molregno = cp.molregno "
            "where md.pref_name is not null "
            "   and cs.canonical_smiles is not null"
            "	and cr.src_id in (1, 8, 9, 12, 13, 36, 41, 42, 53) "
            "and MATCH(md.pref_name) AGAINST (%s IN BOOLEAN MODE) "
            "and cs.canonical_smiles is not null "
            "group by standard_inchi_key, pref_name "
            "order by md.max_phase desc, src_priority ", [query])

        if rows_count > 0:
            rows = cursor.fetchall()
        else:
            return None

    return rows


def check_molecule_synonyms_full_text_search(query):
    with connections['chembl_mysql'].cursor() as cursor:

        rows_count = cursor.execute("select cs.canonical_smiles, cs.standard_inchi, cs.standard_inchi_key, md.max_phase, "
                                    "case "
                                    "  when cr.src_id in (8, 9, 12, 13, 36, 41, 42, 53) then 1 "
                                    "  else 2 "
                                    "end as src_priority, "
                                    "case "
                                    " when md.structure_type = 'NONE' and cp.full_molformula is not null then 'Metal-containing compound' "
                                    " else md.molecule_type end molecule_type, "
                                    "upper(ms.synonyms) "
                                    "from molecule_synonyms ms "
                                    " join compound_records cr on ms.molregno = cr.molregno "
                                    " join molecule_hierarchy mh on cr.molregno = mh.molregno "
                                    " join compound_structures cs on mh.parent_molregno = cs.molregno "
                                    " join molecule_dictionary md on mh.parent_molregno = md.molregno "
                                    " left join compound_properties cp ON mh.parent_molregno = cp.molregno "
                                    "where cr.src_id in (1, 8, 9, 12, 13, 36, 41, 42, 53) "
                                    " and MATCH(ms.synonyms) AGAINST (%s IN BOOLEAN MODE) "
                                    "and cs.canonical_smiles is not null "
                                    "group by standard_inchi_key, synonyms "
                                    "order by md.max_phase desc, src_priority ", [query])

        if rows_count > 0:
            rows = cursor.fetchall()
        else:
            return None

    return rows


def structure_to_pref_name(query):
    structure_to_prefname = '''
    select distinct cs.canonical_smiles, cs.standard_inchi_key,
            md.pref_name,
            case
            when md.structure_type = 'NONE' and cp.full_molformula is not null then 'Metal-containing compound'
            else md.molecule_type end molecule_type,            
            case
            when cr.src_id in (8, 9, 12, 13, 36, 41, 42, 53) then 1
            else 2 end as src_priority,
            md.max_phase
            from compound_structures cs
                left join molecule_dictionary md on cs.molregno = md.molregno
                left join compound_properties cp ON cs.molregno = cp.molregno
                left join compound_records cr ON cs.molregno = cr.molregno
            where md.pref_name is not null
            and cs.standard_inchi_key = %s
            order by md.max_phase desc, src_priority limit 1;
    '''

    with connections['chembl_mysql'].cursor() as cursor:

        rows_count = cursor.execute(structure_to_prefname, [query])

        if rows_count > 0:
            rows = cursor.fetchall()
        else:
            return None

    return rows


def structure_to_synonym(query):
    structure_to_synonym = '''
    select cs.canonical_smiles, cs.standard_inchi_key,
            upper(ms.synonyms), 
            case
            when md.structure_type = 'NONE' and cp.full_molformula is not null then 'Metal-containing compound'
            else md.molecule_type end molecule_type,            
            case
            when cr.src_id in (8, 9, 12, 13, 36, 41, 42, 53) then 1
            else 2
            end as src_priority,
            md.max_phase
            from compound_structures cs
                join molecule_dictionary md on cs.molregno = md.molregno
                left join compound_properties cp ON cs.molregno = cp.molregno
                left join compound_records cr ON cs.molregno = cr.molregno
                left join molecule_synonyms ms ON cs.molregno = ms.molregno
            where cs.standard_inchi_key = %s
            group by cs.canonical_smiles, cs.standard_inchi, cs.standard_inchi_key, md.max_phase, src_priority
            order by md.max_phase desc, src_priority limit 1;
    '''

    with connections['chembl_mysql'].cursor() as cursor:

        rows_count = cursor.execute(structure_to_synonym, [query])

        if rows_count > 0:
            rows = cursor.fetchall()
        else:
            return None

    return rows


class S2NResult:
    def __init__(self, canonical_smiles, standard_inchi_key, name):
        self.name = name
        self.smiles = canonical_smiles
        self.inchikey = standard_inchi_key


class N2SResult:
    def __init__(self, name, smiles, inchikey, truncated_inchikey):
        self.name = name
        self.inchikey = inchikey
        self.smiles = smiles
        self.connectivity = truncated_inchikey


class PAStandardizePreClinicalResult:
    def __init__(self, smiles, inchikey, truncated_inchikey):
        self.inchikey = inchikey
        self.smiles = smiles
        self.connectivity = truncated_inchikey


def parse_n2s_result_rows_for_n2s_query(rows):
    result_list = []

    if rows is None or len(rows) == 0:
        return result_list

    for row in rows:
        name = row[6]
        smiles = row[0]
        inchi_key = row[2]
        truncated_inchikey = get_truncated_inchi_key(inchi_key)

        result_list.append(N2SResult(name, smiles, inchi_key, truncated_inchikey))

    result = {"result": result_list}
    return result


def parse_s2n_result_rows_for_s2n_query(rows):
    result_list = []

    if rows is None or len(rows) == 0:
        return result_list

    for row in rows:
        canonical_smiles = row[0]
        standard_inchi_key = row[1]
        name = row[2]

        result_list.append(S2NResult(canonical_smiles, standard_inchi_key, name))

    result = {"result": result_list}
    return result


def parse_n2s_result_rows_for_pa_standardize_query_as_tuple(rows):
    if rows is None or len(rows) == 0:
        return {}

    for row in rows:
        result = row[2], row[0]
        result_data = {"result": result}
        return result_data
