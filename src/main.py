import pandas as pd
from google.cloud import bigquery
import json

def load_bq(client, project, dataset, table, data, schema):
    '''
    Load table and schema to BigQuery
    '''
    print('Loading: '+dataset+'.'+table)
    
    table_bq = '%s.%s.%s' % (project, dataset, table)
    
    job_config = bigquery.LoadJobConfig( 
        write_disposition="WRITE_TRUNCATE",      
        autodetect=False,
        schema=schema,
        source_format=bigquery.SourceFormat.CSV,
        allow_jagged_rows=True,
        allow_quoted_newlines=True
    )
    
    job = client.load_table_from_dataframe(
        data, table_bq, job_config=job_config
    )


def func(data, context):

    bq_project = 'htan-dcc'
    source_ds = 'combined_assays'
    dest_ds = 'id_provenance'
    dest_table = 'upstream_ids'

    # instantiate BigQuery client
    try:
        client = bigquery.Client()
        print('BigQuery client successfully initialized')
    except:
        print('Failure to initialize BigQuery client')

    bios = client.query("""
        SElECT * FROM `htan-dcc.id_provenance.biospecimen_ids`
        """).result().to_dataframe()

    table_list = client.list_tables(source_ds)

    all_tables = client.query("""
    SELECT table_name,column_name FROM
    `htan-dcc.combined_assays.INFORMATION_SCHEMA.COLUMNS`
    """).result().to_dataframe()

    bq_schema = json.load(open('schema.json'))

    f = pd.DataFrame()

    for t in table_list:

        tn = t.table_id

        if tn not in ['OtherAssay','ExSeqMinimal'] and 'Auxiliary' not in tn and 'Level' not in tn and 'Xenium' not in tn:
            continue

        else:

            print( '' )
            print( ' Processing: ' + str(tn) )
            print( '' )

            cols = list(all_tables[all_tables['table_name'] == tn]['column_name'])

            common_cols = ['HTAN_Data_File_ID', 'Filename', 'entityId', 
                'Component', 'Data_Release', 'CDS_Release', 'HTAN_Center']

            # Add 'HTAN_Parent_Biospecimen_ID' if present in cols
            if 'HTAN_Parent_Biospecimen_ID' in cols:
                common_cols.append('HTAN_Parent_Biospecimen_ID')

            # Add 'HTAN_Parent_Data_File_ID' if present in cols
            if 'HTAN_Parent_Data_File_ID' in cols:
                common_cols.append('HTAN_Parent_Data_File_ID')

            # Join the columns into a string
            qco = ', '.join(common_cols)

            df = client.query("""
                SELECT %s FROM `%s.%s.%s`
                """ % (qco,t.project,t.dataset_id,t.table_id)).result().to_dataframe()

            f = pd.concat([f, df], ignore_index=True)

    # expand comma and semicolon separated parent lists into individual rows
    f = f.assign(HTAN_Parent_Data_File_ID = \
        f.HTAN_Parent_Data_File_ID.str.split("[,;]")).explode('HTAN_Parent_Data_File_ID')

    f = f.assign(HTAN_Parent_Biospecimen_ID = \
        f.HTAN_Parent_Biospecimen_ID.str.split("[,;]")).explode('HTAN_Parent_Biospecimen_ID')

    # trim whitespace
    f = f.applymap(lambda x: x.strip() if isinstance(x, str) else x).drop_duplicates()

    print( '' )
    print( ' Walking parent file ancestry ' )
    print( '' )

    # join on parent IDs. HTAN assays currently have at most 4 levels
    f = (
        f.merge(f[['HTAN_Data_File_ID','HTAN_Parent_Data_File_ID']], 
                left_on='HTAN_Parent_Data_File_ID',
                right_on='HTAN_Data_File_ID',
                how='left')
        .drop_duplicates()
        .rename(columns={'HTAN_Parent_Data_File_ID_y':'gParent_File_ID',
                        'HTAN_Parent_Data_File_ID_x':'HTAN_Parent_Data_File_ID',
                        'HTAN_Data_File_ID_x':'HTAN_Data_File_ID'})
        .drop(columns='HTAN_Data_File_ID_y')
        .merge(f[['HTAN_Data_File_ID','HTAN_Parent_Data_File_ID']], 
                left_on='gParent_File_ID',
                right_on='HTAN_Data_File_ID',
                how='left')
        .drop_duplicates()
        .rename(columns={'HTAN_Data_File_ID_x':'HTAN_Data_File_ID',
                        'HTAN_Parent_Data_File_ID_y':'ggParent_File_ID',
                        'HTAN_Parent_Data_File_ID_x':'HTAN_Parent_Data_File_ID'})
        .drop(columns='HTAN_Data_File_ID_y')
    )

    # coalesce 'highest' level parent data file
    f['coalesce'] = f[['ggParent_File_ID','gParent_File_ID', \
        'HTAN_Parent_Data_File_ID']].bfill(axis=1).iloc[:,0]

    f = (
        f.merge(f[['HTAN_Data_File_ID', 'HTAN_Parent_Biospecimen_ID']], 
                left_on='coalesce',
                right_on='HTAN_Data_File_ID',
                how='left')
        .drop_duplicates()
        .assign(HTAN_Parent_Biospecimen_ID_x=lambda x: x[
            'HTAN_Parent_Biospecimen_ID_x'].fillna(x['HTAN_Parent_Biospecimen_ID_y']))
        .rename(columns={'HTAN_Data_File_ID_x':'HTAN_Data_File_ID',
                        'HTAN_Parent_Biospecimen_ID_x':'HTAN_Assayed_Biospecimen_ID'})
        .drop(columns=['HTAN_Data_File_ID_y', 'gParent_File_ID', 'coalesce', 
            'HTAN_Data_File_ID_y', 'HTAN_Parent_Biospecimen_ID_y', 'ggParent_File_ID'])
    )

    # join file table with biospecimen walker table
    id_prov = f.merge(bios, 
        on='HTAN_Assayed_Biospecimen_ID', how='left').drop_duplicates()

    load_bq(client, bq_project, dest_ds, dest_table, id_prov, bq_schema)

    print( '' )
    print( ' Done! ')
    print( '' )

