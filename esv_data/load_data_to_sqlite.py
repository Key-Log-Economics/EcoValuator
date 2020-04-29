"""Loads data from CSV files to SQLite database
Expects CSV files with columns in the given order:

esv_estimates.csv
    lulc_source
    lulc_value
    service_id
    estimate_min
    estimate_max
    estimate_avg

lulc_legend.csv
    source
    value
    description

service_names.csv
    service_id
    service_name

run this file within the directory of the csv files to be loaded into a database file
THIS SCRIPT WILL OVERWRITE EXISTING DATABASE FILE IN THE DIRECTORY
"""
import os
import sqlite3
import csv

DATA_DIR = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

DB_FILE = os.path.join(DATA_DIR, 'ESV_data.sqlite')

if os.path.exists(DB_FILE):
    os.remove(DB_FILE)

SQL_INFO = [{'table_name':'esv_estimates',
             'columns': [('lulc_source','TEXT','NOT NULL'),
                        ('lulc_value','INTEGER','NOT NULL'), 
                        ('service_id','INTEGER','NOT NULL'),
                        ('estimate_min','REAL','NOT NULL'),
                        ('estimate_max','REAL','NOT NULL'),
                        ('estimate_avg','REAL','NOT NULL')]
            },
            {'table_name':'lulc_legend',
             'columns':[('source', 'TEXT','NOT NULL'),
                        ('value','INTEGER','NOT NULL'),
                        ('description','TEXT','NOT NULL')]
            },
            {'table_name':'service_names',
            'columns':[('service_id','INTEGER','NOT NULL'),
                       ('service_name','TEXT','NOT NULL')]
            }]


for table in SQL_INFO:
    table_name = table['table_name']
    print(f'processing table {table_name} ')
    columns = table['columns']
    column_names = [c[0] for c in columns]
    columns_sql = [' '.join(col) for col in columns]

    qs = '?,' * len(columns)
    qs = qs[:-1]

    create_sql = f"""CREATE TABLE IF NOT EXISTS {table_name}({','.join(columns_sql)})"""
    insert_sql = f"""INSERT INTO {table_name}({','.join(column_names)}) VALUES ({qs})"""

    csv_file = os.path.join(DATA_DIR, table_name + '.csv')
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        rows = list(reader)
        columns = rows[0]
        data = rows[1:]

    with sqlite3.connect(DB_FILE) as conn:
        cursor = conn.cursor()
        cursor.execute(create_sql)
        cursor.executemany(insert_sql, data)

