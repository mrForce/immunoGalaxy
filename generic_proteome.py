import DB
import argparse

desc = """Annoyingly, when I created the tidePipeline projects for Human and Mouse, I named the peptide lists and NetMHC runs "Human*" and "Mouse*". Now, I want to genericize those names. I think that the easiest way to do this will be to simply change the names in the existing sqlite database.

Set PYTHONPATH so that we can import DB, give the path to database.db, and give the prefix to change, and what to change it to. For example, pass Mouse and Proteome so "Mouse8Mers" is changed to "Proteome8Mers"
"""

parser = argparse.ArgumentParser(description=desc)
parser.add_argument('db_path', help='path to database.db')
parser.add_argument('prefix', help='Prefix to match and change in PeptideList and NetMHC rows')
parser.add_argument('change_to', help='What to change the prefix to')

args = parser.parse_args()
session = DB.init_session(args.db_path)
for peptidelist_row in session.query(DB.PeptideList).all():
    assert(peptidelist_row.peptideListName.startswith(args.prefix))

for netmhc_row in session.query(DB.NetMHC).all():
    assert(netmhc_row.Name.startswith(args.prefix))

for peptidelist_row in session.query(DB.PeptideList).all():
    peptidelist_row.peptideListName = peptidelist_row.peptideListName.replace(args.prefix, args.change_to, 1)

for netmhc_row in session.query(DB.NetMHC).all():
    netmhc_row.Name = netmhc_row.Name.replace(args.prefix, args.change_to, 1)

session.commit()

