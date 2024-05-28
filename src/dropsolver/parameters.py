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
        for field in ['dimension', 'html', 'display', 'comments']:
            if not parameter[field]:
                del parameter[field]
        parameters.append(parameter)
    return parameters

def to_defaults_string():
    defaults = []
    for parameter in parameters():
        if 'default_value' not in parameter:
            continue
        if 'SI_multiplier' in parameter:
            parameter['default_value'] *= parameter['SI_multiplier']
        defaults.append("{}={}".format(parameter['parameter'], parameter['default_value']))
    return ', '.join(defaults)
