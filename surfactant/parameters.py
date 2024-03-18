import csv
import re

def parameters():
    with open('parameters.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        lines = list(reader)
    header = lines.pop(0)
    header = [re.sub(' ', '_', field) for field in header]
    parameters = []
    for fields in lines:
        parameter = dict(zip(header, fields))
        for field in ['SI_multiplier', 'default_value', 'min', 'max']:
            if parameter[field]:
                parameter[field] = float(parameter[field])
            else:
                del parameter[field]
        for field in ['dimension', 'comments']:
            if not parameter[field]:
                del parameter[field]
        parameters.append(parameter)
    return parameters
