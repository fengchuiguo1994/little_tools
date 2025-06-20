import csv
import logging
import multiprocessing as mp

from contextlib import contextmanager
from pathlib import Path
from typing import Dict
from mysql.connector import MySQLConnection
from mysql.connector.cursor import MySQLCursor

import mysql.connector
import random

logging.basicConfig(level=logging.DEBUG)

config = {
    'host': "localhost",
    'username': "root",
    'password': "kkltmax8",
    'database': "test",
    'autocommit': False
}

@contextmanager
def mysql_connect(connection_config: Dict) -> MySQLConnection:
    conn: MySQLConnection = mysql.connector.connect(**connection_config)
    # logging.info('Connected to database \'%s\' at \'%s\'', connection_config['database'], connection_config['host'])
    yield conn

    logging.info('Closing connection to server \'%s\'', connection_config['host'])
    conn.close()

@contextmanager
def mysql_cursor(conn: MySQLConnection) -> MySQLCursor:
    cur: MySQLCursor = conn.cursor()
    yield cur
    # logging.info('Closing cursor')
    cur.close()


insert_sql = """INSERT INTO parking_citations (col1, col2, col3) VALUES ( %s, %s, %s ) """

def worker(rows):
    with mysql_connect(config) as conn:
        with mysql_cursor(conn) as cursor:
            logging.info('Processing batch')
            cursor.executemany(insert_sql, rows)
            conn.commit()

def get_chunks(batch_size: int, source_file: Path):
    batch_data = []
    batch_count = 0
    for i in range(10000000):    
        batch_data.append([f"value_{batch_count}", random.randint(1, 100), random.uniform(1.0, 100.0)])
        batch_count += 1
        if batch_count % batch_size == 0:
            yield batch_data
            batch_data = []
    if batch_data:
        yield batch_data

def main():
    batch_size = 5000
    source_file = Path('data/parking-citations.csv')
    chunk_gen = get_chunks(batch_size, source_file)
    # pool = mp.Pool(mp.cpu_count()-1)
    pool = mp.Pool(8)
    print(mp.cpu_count())
    results = pool.imap(worker, chunk_gen)
    pool.close()
    pool.join()

if __name__ == '__main__':
    host = "localhost"
    user = "root"
    password = "kkltmax8"
    spe = "assembly"
    database = "test"

    db = mysql.connector.connect(
        host = host,
        user = user,
        password = password,
        database = database,
        connect_timeout=60
    )

    cursor = db.cursor()

    # 创建表
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS parking_citations (
            id INT AUTO_INCREMENT PRIMARY KEY,
            col1 VARCHAR(50),
            col2 INT,
            col3 FLOAT
        )
    """)
    db.close()

    main()