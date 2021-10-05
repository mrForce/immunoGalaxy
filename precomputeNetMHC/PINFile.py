import csv
import re
def parse_peptide(peptide, peptide_regex, ptm_removal_regex):
    match = peptide_regex.match(peptide)
    if match and match.group('peptide'):
        matched_peptide = match.group('peptide')
        return ptm_removal_regex.sub('', matched_peptide)
    else:
        return None

class PINFile:
    def __init__(self, path):
        self.path = path
        self.peptides = set()
        self.peptide_regex = re.compile('^[A-Z\-]\.(?P<peptide>.*)\.[A-Z\-]$')
        self.ptm_removal_regex = re.compile('\[[^\]]*\]')
        with open(path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t', restkey='Proteins')
            self.fieldnames = list(reader.fieldnames)
            next(reader)
            for row in reader:
                peptide = parse_peptide(row['Peptide'], self.peptide_regex, self.ptm_removal_regex)
                assert(peptide)
                self.peptides.add(peptide)
    def addPin(self, additionalPath):
        rows = []
        with open(additionalPath, 'r') as f:
            reader = csv.DictReader(f, self.fieldnames, delimiter='\t', restkey='Proteins')
            next(reader)
            for row in reader:
                rows.append(row)
            with open(self.path, 'w') as g:
                writer = csv.DictWriter(g, self.fieldnames, delimiter='\t')
                for row in rows:
                    writer.writerow(row)
                    
            
    def addScores(self, scoreDict, columnHeader, columnDirection):
        #inserts the scores
        rows = []
        with open(self.path, 'r') as f:
            reader = csv.DictReader(f, self.fieldnames, delimiter='\t', restkey='Proteins')
            next(reader)
            direction = dict(next(reader))
            for row in reader:
                rows.append(row)
        self.fieldnames.insert(6, columnHeader)
        direction[columnHeader] = columnDirection
        with open(self.path, 'w') as f:
            writer = csv.DictWriter(f, self.fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerow(direction)
            for row in rows:
                row_copy = dict(row)
                peptide = parse_peptide(row['Peptide'], self.peptide_regex, self.ptm_removal_regex)
                row_copy[columnHeader] = str(scoreDict[peptide])
                writer.writerow(row_copy)
