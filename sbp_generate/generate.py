#!/usr/bin/env python

import yaml
import jinja2
import re

with open("sbp.yaml", 'r') as f:
  ds = yaml.load(f)

LATEX_SUBS = (
    (re.compile(r'\\'), r'\\textbackslash'),
    (re.compile(r'([{}_#%&$])'), r'\\\1'),
    (re.compile(r'~'), r'\~{}'),
    (re.compile(r'\^'), r'\^{}'),
    (re.compile(r'_'), r'_'),
    (re.compile(r'"'), r"''"),
    (re.compile(r'\.\.\.+'), r'\\ldots'),
)

def escape_tex(value):
    newval = value
    for pattern, replacement in LATEX_SUBS:
        newval = pattern.sub(replacement, newval)
    return newval

def commentify(value):
  return '\n'.join([' * ' + l for l in value.split('\n')[:-1]])

acronyms = ['GPS', 'ECEF', 'LLH', 'NED']
def classnameify(s):
  return ''.join(w if w in acronyms else w.title() for w in s.split('_'))

jenv = jinja2.Environment(
    block_start_string = '((*',
    block_end_string = '*))',
    variable_start_string = '(((',
    variable_end_string = ')))',
    comment_start_string = '((=',
    comment_end_string = '=))',
    loader=jinja2.FileSystemLoader("./")
)
jenv.filters['escape_tex'] = escape_tex
jenv.filters['commentify'] = commentify
jenv.filters['classnameify'] = classnameify

sizes = {
    'u8': 1,
    'u16': 2,
    'u32': 4,
    'u64': 8,
    's8': 1,
    's16': 2,
    's32': 4,
    's64': 8,
    'float': 4,
    'double': 8,
}

pystruct_code = {
    'u8': 'B',
    'u16': 'H',
    'u32': 'I',
    'u64': 'Q',
    's8': 'b',
    's16': 'h',
    's32': 'i',
    's64': 'q',
    'float': 'f',
    'double': 'd',
}


def pystruct_format(fields):
  return '<' + ''.join(pystruct_code[f['type']] for f in fields)
jenv.filters['pystruct'] = pystruct_format

def rejig_values(values):
    new_values = []
    for v in values:
        value, desc = v.iteritems().next()
        new_values.append({
            'value': value,
            'desc': desc
        })
    return new_values

def rejig_bitfields(bfs):
    new_bfs = []
    n_with_values = 0
    for bf in bfs:
        rng, info = bf.iteritems().next()
        if 'values' in info:
            n_with_values += 1
            info['vals'] = rejig_values(info['values'])
            del info['values']
        rng = map(int, str(rng).split('-'))
        if len(rng) == 1:
            lsb = rng[0]
            bf_len = 1
        else:
            lsb, msb = rng
            bf_len = msb + 1 - lsb
        rng = ':'.join(map(str, rng))
        new_bfs.append(dict(info, **{
            'lsb': lsb,
            'range': rng,
            'len': bf_len,
        }))
    return new_bfs, n_with_values

msgs = [[dict({'name': k}, **v) for k, v in d.iteritems()][0] for d in ds]
for m in msgs:
    fields = []
    offset = 0
    max_type_len = 0
    max_name_len = 0
    for f in m['fields']:
        name, info = f.iteritems().next()
        if not 'units' in info:
            info.update({'units': ''})
        info['n_with_values'] = 0
        if 'fields' in info:
            info['fields'], info['n_with_values'] = rejig_bitfields(info['fields'])
        fields.append(dict(info, **{
            'name': name,
            'offset': offset,
            'size': sizes[info['type']],
        }))
        offset += sizes[info['type']]
        max_type_len = max(max_type_len, len(info['type']))
        max_name_len = max(max_name_len, len(name))
    m['max_type_len'] = max_type_len
    m['max_name_len'] = max_name_len
    m['fields'] = fields
    m['size'] = offset

latex_template = jenv.get_template('message_descs.tex')
with open("sbp_out.tex", 'w') as f:
    f.write(latex_template.render(msgs=msgs))

c_template = jenv.get_template('sbp_messages_template.h')
with open("../include/libswiftnav/sbp_messages.h", 'w') as f:
    f.write(c_template.render(msgs=msgs))

py_template = jenv.get_template('sbp_messages_template.py')
with open("sbp_messages.py", 'w') as f:
    f.write(py_template.render(msgs=msgs))

import subprocess

subprocess.call(["pdflatex" , "sbp_out.tex"])
subprocess.call(["mv" , "sbp_out.pdf", "../docs/sbp.pdf"])

