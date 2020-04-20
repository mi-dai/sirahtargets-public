import sqlite3
from sqlite3 import Error
import pandas as pd

def create_connection(db_file):
    """ create a database connection to a SQLite database """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
    except Error as e:
        print(e)
    finally:
        if conn:
            conn.close()
            
# create_connection('db/glade_v23.db')
df = pd.read_csv('catalogues/GLADE/GLADE_2.3.txt',header=None,sep='\s+')
con = sqlite3.connect("db/glade_v23.db")

# drop data into database
df.to_sql("glade_v23", con, if_exists='replace')

con.close()