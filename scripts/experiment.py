#!/usr/bin/env python
# coding: utf-8

import json

# Export experiment information to json file.
diction = {'sample': 'ds193',
           'names': ['4K 100nJ', '4K 10nJ', '4K 350nJ', '4K 1nJ', '50K 100nJ', '50K 10nJ', '50K 425nJ', '50K 1nJ', 
                     '200K 100nJ', '200K 10nJ', '200K 390nJ', '200K 1nJ', '300K 390nJ'],
           'fluences': [100, 10, 350, 1, 100, 10, 425, 1, 100, 10, 390, 1, 390],
           'temperatures': [4, 4, 4, 4, 50, 50, 50, 50, 200, 200, 200, 200, 300],
           'integration': [5, 5, 2, 10, 5, 15, 5, 10, 5, 15, 5, 10, 5]}

with open('experiment.json', 'w', encoding='utf-8') as f:
    json.dump(diction, f, ensure_ascii=False, indent=4)


# Test import of file into string and list.
with open('experiment.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    search = data['sample']
    names = data['names']
    fluences = data['fluences']
    temperatures = data['temperatures']

print(search)
print(names)
